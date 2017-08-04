
#include <unistd.h>
#include <sys/time.h>

#include "db.hpp"
#include "parse.hpp"


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


PairDB::PairDB(const char * path) : m_path(path), m_worker_thread(0), m_head(0), m_tail(0), m_num_written(0)
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

PairDB::~PairDB()
{
    close();
}

void
PairDB::close()
{
    if (NULL != m_handle) {
        sqlite3_close(m_handle);
        m_handle = NULL;
    }
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
                           "VALUES (1, %d, 1, %Q, %d, %d, 1, %Q)",
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
    if (0 != m_worker_thread)
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
                                   "VALUES (1, %d, 1, %Q, %d, %d, 1, %Q)",
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

void
PairDB::store_counters(Counters * c, int cotrans_min_length)
{
    SQL_LOCAL_VARS(m_handle);
    SQL_EXEC("CREATE TABLE IF NOT EXISTS counter (run_key TEXT, dict_index INT, count_key TEXT, count INT)");
    SQL_EXEC("CREATE INDEX IF NOT EXISTS counter_key_idx ON counter (run_key)");
    SQL_EXEC("BEGIN");
    char key[1024] = { 0 };
    for (int mask = MASK_TREATED; mask <= MASK_UNTREATED; ++mask) {
        for (int L = cotrans_min_length; L < c->n; ++L) {
            for (int site = 0; site <= L; ++site) {
                snprintf(key, 1023, "1:%s:%d:%d", mask_chars(mask), site, L);
                zSql = sqlite3_mprintf("INSERT INTO counter VALUES ('spats', 1, %Q, %d)", key, c->site_count(mask, L, site));
                SQL_EXEC(zSql);
            }
        }
    }
    SQL_EXEC("COMMIT");
}

void
PairDB::store_run()
{
    SQL_LOCAL_VARS(m_handle);
    SQL_EXEC("CREATE TABLE IF NOT EXISTS run_data (param_key TEXT, param_val TEXT)");
    SQL_EXEC("DELETE FROM run_data");
    SQL_EXEC("INSERT INTO run_data VALUES ('algorithm', 'native')");
    SQL_EXEC("INSERT INTO run_data VALUES ('cotrans', 'True')");
}

void
PairDB::store_targets(Targets * targets)
{
    SQL_LOCAL_VARS(m_handle);
    SQL_EXEC("CREATE TABLE IF NOT EXISTS target (name TEXT, seq TEXT)");
    SQL_EXEC("DELETE FROM target");
    for (int i = 0; i < targets->size(); ++i) {
        Target * t = targets->target(i);
        zSql = sqlite3_mprintf("INSERT INTO target VALUES (%Q, %Q)", t->name().c_str(), t->seq().c_str());
        SQL_EXEC(zSql);
    }
}

struct ParseInfo
{
    PairDB * db;
    int appx_total;
    int parsed;
    int left_to_sample;
    sqlite3_stmt * stmt;
};

ParseInfo * g_parsing_info = NULL;

bool
parse_db_handler(Fragment * r1, Fragment * r2, const char * handle)
{
    ParseInfo * pi = g_parsing_info;
    ++pi->parsed;

    long r = random();
    if (pi->parsed >= pi->appx_total) {
        if (pi->left_to_sample == 0)
            return false;
    }
    else if (r % (pi->appx_total - pi->parsed) >= pi->left_to_sample) {
        return true;
    }

    std::string str(handle);
    str += r1->string();
    int sqlErr = sqlite3_bind_text(pi->stmt, 1, str.c_str(), (int)str.length(), SQLITE_TRANSIENT); // SQLITE_STATIC?
    if (sqlErr != SQLITE_OK)
        return false;

    str = r2->string();
    sqlErr = sqlite3_bind_text(pi->stmt, 2, str.c_str(), (int)str.length(), SQLITE_TRANSIENT); // SQLITE_STATIC?
    if (sqlErr != SQLITE_OK)
        return false;

    sqlErr = sqlite3_step(pi->stmt);
    if (sqlErr != SQLITE_DONE)
        return false;

    sqlite3_reset(pi->stmt);

    if (0 == --pi->left_to_sample)
        return false;

    return true;
}

void
PairDB::parse_and_sample(const char * r1_path, const char * r2_path, int sample_size)
{
    SQL_LOCAL_VARS(m_handle);
    SQL_EXEC("CREATE TABLE IF NOT EXISTS pair (r1 TEXT, r2 TEXT, identifier TEXT)");
    SQL_EXEC("DELETE FROM pair");
    SQL_EXEC("BEGIN");

    ATS_DEBUG("SS: %d\n", sample_size);

    portable_srandomdev();

    ParseInfo pi;
    pi.db = this;
    pi.appx_total = appx_number_of_fastq_pairs(r1_path);
    pi.parsed = 0;
    pi.left_to_sample = sample_size;

    sqlErr = sqlite3_prepare_v2(m_handle, "INSERT INTO pair (identifier, r1, r2) VALUES ('', ?, ?)", -1, &pi.stmt, NULL);
    ATS_THROW_IF(sqlErr != SQLITE_OK || NULL == pi.stmt);

    ATS_ASSERT(NULL == g_parsing_info);
    g_parsing_info = &pi;
    fastq_parse_handler(r1_path, r2_path, &parse_db_handler);
    g_parsing_info = NULL;

    SQL_EXEC("COMMIT");
}
