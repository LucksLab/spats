
#include <unistd.h>
#include <sys/time.h>

#include "db.hpp"

typedef struct
{
    int (*callback)(void*,int,char**,char**);
    void * param;
} _SM_CPP_CallbackInfo_t;

/* for internal use -- to make sure exceptions are always caught in callbacks */
int
_SM_cppWrappingCallback(void * param, int numColumns, char** colValues, char** colNames);
int
_SM_cppWrappingFailOKCallback(void * param, int numColumns, char** colValues, char** colNames);


#define SQL_LOCAL_VARS(theHandle) ATS_THROW_IF(NULL == (theHandle)); sqlite3 * activeDbHandle = (theHandle); char * zSql = NULL; int sqlErr; char * errorMessage = NULL
#define SQL_CHANGE_HANDLE(theDbHandle) ATS_THROW_IF(NULL == (theDbHandle)); activeDbHandle = (theDbHandle);

#define DO_SQL_EXEC(query, errval, theCallback, theParam, errorLogger, debugLogger,errorHandler,callbackWrapper) do { \
        debugLogger("%s", query);                                       \
        {                                                               \
            _SM_CPP_CallbackInfo_t callbackinfo;                        \
            callbackinfo.callback = theCallback;                        \
            callbackinfo.param = theParam;                              \
            sqlErr = sqlite3_exec(activeDbHandle, query, &callbackWrapper, &callbackinfo, &errorMessage); \
        }                                                               \
        if (sqlErr != SQLITE_OK) {                                      \
            errorLogger("sql error: %s at %s:%d\nquery: %s", errorMessage, __FILE__, __LINE__, query); \
            ATS_DEBUG("  at: ");                                        \
            if (errorMessage != NULL) sqlite3_free(errorMessage);       \
            errorMessage = NULL;                                        \
            errorHandler;                                               \
        }                                                               \
        if (errorMessage != NULL) sqlite3_free(errorMessage);           \
        errorMessage = NULL;                                            \
        if (zSql != NULL) sqlite3_free(zSql);                           \
        zSql = NULL;                                                    \
    } while (0)


#define SM_SQLMACRO_IGNORE(...)

#ifndef SM_SQL_EXEC_DEBUG_LOGGER
# if AU_SM_DEBUG_TRACE_SQL
#  define SM_SQL_EXEC_DEBUG_LOGGER ATS_DEBUG
# else
#  define SM_SQL_EXEC_DEBUG_LOGGER SM_SQLMACRO_IGNORE
# endif
#endif

#ifndef SM_SQL_EXEC_ERROR_LOGGER
# define SM_SQL_EXEC_ERROR_LOGGER ATS_WARN
#endif

#define SQL_EXEC_LOG(query,errval,callback,param,errorHandler,callbackWrapper) DO_SQL_EXEC(query,errval,callback,param,SM_SQL_EXEC_ERROR_LOGGER,SM_SQL_EXEC_DEBUG_LOGGER,errorHandler,callbackWrapper)

#define SQL_EXEC(query) SQL_EXEC_LOG(query,NULL,NULL,NULL,throw std::exception(),_SM_cppWrappingCallback)
#define SQL_EXEC_CB(query,callback,param) SQL_EXEC_LOG(query,errval,callback,param,throw std::exception(),_SM_cppWrappingCallback)
#define SQL_EXEC_FAIL_OK(query) SQL_EXEC_LOG(query,NULL,NULL,NULL,;,_SM_cppWrappingCallback)
#define SQL_EXEC_CB_FAIL_OK(query,callback,param) SQL_EXEC_LOG(query,errval,callback,param,;,_SM_cppWrappingFailOKCallback)



int
_SM_cppWrappingCallback(void * param, int numColumns, char** colValues, char** colNames)
{
    _SM_CPP_CallbackInfo_t * callbackInfo = (_SM_CPP_CallbackInfo_t *)param;
    ATS_RETVAL_IF(NULL == callbackInfo, -1);
    try
    {
        return (callbackInfo->callback)(callbackInfo->param, numColumns, colValues, colNames);
    }
    catch (...)
    {
        ATS_WARN("Caught exception in sqlite callback (bp at sqliteHelpers.cpp::_SM_cppWrappingCallback).");
        return -1;
    }
}


int
_SM_cppWrappingFailOKCallback(void * param, int numColumns, char** colValues, char** colNames)
{
    _SM_CPP_CallbackInfo_t * callbackInfo = (_SM_CPP_CallbackInfo_t *)param;
    ATS_RETVAL_IF(NULL == callbackInfo, -1);
    try
    {
        return (callbackInfo->callback)(callbackInfo->param, numColumns, colValues, colNames);
    }
    catch (...)
    {
        ATS_WARN("Caught exception in sqlite callback (bp at sqliteHelpers.cpp::_SM_cppWrappingCFailOKCallback).");
        ATS_ASSERT_NOT_REACHED();
        return 0; 
    }
}


PairDB::PairDB(const char * path) : m_path(path), m_worker_thread(NULL), m_head(0), m_tail(0), m_num_written(0)
{
    int sqliteErr = SQLITE_OK;
    sqliteErr = sqlite3_open_v2(path,
                                &m_handle,
                                (SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE),
                                NULL);
    if (sqliteErr != SQLITE_OK) {
        m_handle = NULL;
        throw std::exception();
    }
    memset(m_cases, 0, sizeof(m_cases));
}

int
test_dump_CB(void * param, int numColumns, char** colValues, char** colNames)
{
    printf("%s / %s\n", colValues[0], colValues[1]);
    return 0; // ret -1 on error
}

void
PairDB::test()
{
    SQL_LOCAL_VARS(m_handle);
    //SQL_EXEC("CREATE TABLE IF NOT EXISTS tester (value INT)");
    //SQL_EXEC("INSERT INTO tester VALUES ('18')");
    //zSql = sqlite3_mprintf("INSERT INTO _bf_SMVersion VALUES(%d, %d)", provider->providerType(), targetVersion);
    //SQL_EXEC(zSql);
    //SQL_EXEC_CB("SELECT version FROM version LIMIT 1", &SM_selectIntCallback, &version);
    SQL_EXEC_CB("SELECT r1, r2 FROM pair LIMIT 5", &test_dump_CB, NULL);
}

int
run_case_CB(void * param, int numColumns, char** colValues, char** colNames)
{
    Spats * spats = (Spats *)param;
    Case c(atoi(colValues[0]), colValues[1], colValues[2]);
    spats->run_case(&c);
    return 0; // ret -1 on error
}

void
PairDB::run_cases(Spats * spats)
{
    SQL_LOCAL_VARS(m_handle);
    SQL_EXEC_CB("SELECT rowid, r1, r2 FROM pair", &run_case_CB, (void *)spats);
}

void
PairDB::write_result(Case * c)
{
    SQL_LOCAL_VARS(m_handle);
    zSql = sqlite3_mprintf("INSERT INTO result (set_id, pair_id, target, mask, site, end, multiplicity, failure) "
                           "VALUES (1, ?, 1, ?, ?, ?, 1, ?)",
                           c->pair_id, mask_chars(c->mask), c->site, c->L, "");
    SQL_EXEC(zSql);
}

void *
db_worker_fn(void * arg)
{
    PairDB * db = (PairDB *)arg;
    db->worker_fn();
    return NULL;
}

void
PairDB::start_worker()
{
    if (NULL != m_worker_thread)
        return;
    pthread_mutex_init(&m_mutex, NULL);
    m_working = true;
    pthread_create(&m_worker_thread, NULL, &db_worker_fn, (void *)this);
}

void
PairDB::worker_fn()
{
    SQL_LOCAL_VARS(m_handle);
    SQL_EXEC("BEGIN");
    while (m_working) {

        while (m_head != m_tail  &&  m_cases[m_head].pair_id > 0) {
            Case * c = &m_cases[m_head];
            zSql = sqlite3_mprintf("INSERT INTO result (set_id, pair_id, target, mask, site, end, multiplicity, failure) "
                                   "VALUES (1, ?, 1, ?, ?, ?, 1, ?)",
                                   c->pair_id, mask_chars(c->mask), c->site, c->L, "");
            SQL_EXEC(zSql);
            ++m_num_written;
            c->pair_id = 0;
            m_head = (m_head + 1) % CASE_QUEUE_LEN;
        }

        usleep(100); // this is highly tune-able based on how quickly work keeps up...

    }
    SQL_EXEC("COMMIT");
}

void
PairDB::submit_result(Case * c)
{
    while ((m_tail + 1) % CASE_QUEUE_LEN == m_head) {
        printf("x");
        usleep(100);
    }
    int tail = 0;
    pthread_mutex_lock(&m_mutex);
    {
        tail = m_tail;
        m_tail = (m_tail + 1) % CASE_QUEUE_LEN;
    }
    pthread_mutex_unlock(&m_mutex);
    Case * c2 = &m_cases[tail];
    c2->mask = c->mask;
    c2->L = c->L;
    c2->site = c->site;
    c2->pair_id = c->pair_id;
}


void
PairDB::commit_results()
{
    m_working = false;
    pthread_join(m_worker_thread, NULL);
}
