
#include <Python.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <string.h>

char * g_sequence;
std::map<std::string, std::vector<int> > g_index;
int g_word_len = 8;

static PyObject *
spats_tindex(PyObject *self, PyObject *args)
{
    const char *command;
    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    int input_len = (int)strlen(command);
    g_sequence = (char *)malloc(input_len + 1);
    strncpy(g_sequence, command, input_len);
    g_sequence[input_len] = 0;
    char * word = (char *)malloc((g_word_len + 1) * sizeof(char));
    for (int i = 0; i < input_len - g_word_len; ++i) {
        strncpy(word, &command[i], g_word_len);
        word[g_word_len] = 0;
        g_index[word].push_back(i);
    }
    free(word);
    return Py_BuildValue("i", 0);
}

static PyObject *
spats_find(PyObject *self, PyObject *args)
{
    const char *command;
    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    return Py_BuildValue("i", g_index[command].at(0));
}

typedef struct
{
    int query_index;
    int match_len;
    int target_index;
} Candidate;


const int nuc_A = (1);
const int nuc_C = (1<<1);
const int nuc_G = (1<<2);
const int nuc_T = (1<<3);

inline int
char_to_mask(char ch)
{
    switch (ch) {
    case 'A' : return nuc_A;
    case 'C' : return nuc_C;
    case 'G' : return nuc_G;
    case 'T' : return nuc_T;
    case 'U' : return nuc_T;
    case 'Y' : return (nuc_C | nuc_T);
    case 'R' : return (nuc_G | nuc_A);
    case 'S' : return (nuc_C | nuc_G);
    case 'W' : return (nuc_A | nuc_T);
    case 'K' : return (nuc_G | nuc_T);
    case 'M' : return (nuc_A | nuc_C);
    case 'B' : return (nuc_C | nuc_G | nuc_T);
    case 'D' : return (nuc_A | nuc_G | nuc_T);
    case 'H' : return (nuc_A | nuc_C | nuc_T);
    case 'V' : return (nuc_A | nuc_C | nuc_G);
    case 'N' : return (nuc_A | nuc_C | nuc_G | nuc_T);
    default:
        return 0;
    }
}

void
longest_match(const char * s1, int i1, int len1, const char * s2, int i2, int len2, int * left_out, int * right_out)
{
    int left1 = i1;
    int left2 = i2;
    int lmax = std::min(left1, left2);
    int right1 = left1 + len1;
    int right2 = left2 + len2;
    int s1len = (int)strlen(s1);
    int s2len = (int)strlen(s2);
    int rmax = std::min(s1len - right1, s2len - right2);

    int left = 1;
    int right = 0;
    char c1 = 0;
    char c2 = 0;
    int m1 = 0;
    int m2 = 0;
    while (1) {
        if (left > lmax)
            break;
        c1 = s1[left1 - left];
        c2 = s2[left2 - left];
        if (c1 != c2) {
            m1 = char_to_mask(c1);
            m2 = char_to_mask(c2);
            if (0 == (m1 & m2))
                break;
        }
        ++left;
    }

    while (1) {
        if (right >= rmax)
            break;
        c1 = s1[right1 + right];
        c2 = s2[right2 + right];
        if (c1 != c2) {
            m1 = char_to_mask(c1);
            m2 = char_to_mask(c2);
            if (0 == (m1 & m2))
                break;
        }
        ++right;
    }
    *left_out = left - 1;
    *right_out = right;
    return;
}

static PyObject *
spats_find_partial(PyObject *self, PyObject *args)
{
    const char * query;
    int min_len;
    if (!PyArg_ParseTuple(args, "si", &query, &min_len))
        return NULL;
    if (min_len < g_word_len) {
        Py_RETURN_NONE;
    }
    int check_every = min_len - g_word_len;
    if (check_every <= 0)
        check_every = 1;
    int query_len = (int)strlen(query);
    Candidate candidate = { 0 };
    char * site_key = (char *)malloc((g_word_len + 1) * sizeof(char));
    int left, right, index, total_len;

    // NOTE: it's important to check all sites, and all hits -- to find the longest match.
    for (int site = 0; site < query_len; site += check_every) {
        if (site > query_len - check_every)
            site = query_len - site;
        strncpy(site_key, &query[site], g_word_len);
        site_key[g_word_len] = 0;
        std::vector<int>::iterator it;
        for (it = g_index[site_key].begin(); it != g_index[site_key].end(); ++it) {
            index = *it;
            longest_match(query, site, g_word_len, g_sequence, index, g_word_len, &left, &right);
            total_len = left + right + g_word_len;
            if (total_len >= min_len) {
                if (total_len > candidate.match_len) {
                    candidate.query_index = site - left;
                    candidate.match_len = total_len;
                    candidate.target_index = index - left;
                }
                if (total_len >= query_len)
                    goto find_partial_done;
            }
        }
    }
find_partial_done:
    free(site_key);
    return Py_BuildValue("iii", candidate.query_index, candidate.match_len, candidate.target_index);
/*
        for site in check_sites:
            #print "CS: {}, {}".format(site, site_key)
            for index in self._index.get(site_key, []):
                #print "GOT: " + str(index)
                left, right = longest_match(query, (site, word_len), self.seq, (index, word_len))
                total_len = left + right + word_len
                #print "extends: <--{}, -->{} / {} ({})".format(left, right, total_len, min_len)
                if total_len >= min_len:
                    if total_len >= query_len:
                        # we can return immediately if we've got a full match
                        return site - left, total_len, index - left
                    elif not candidate[1] or total_len > candidate[1]:
                        # ...otherwise, keep it if it's the best match so far
                        candidate = (site - left, total_len, index - left)
                        #print "C: {}".format(candidate)
        return candidate
*/
}

static PyMethodDef SpatsMethods[] = {
    {"tindex",  spats_tindex, METH_VARARGS, "Index test."},
    {"find",  spats_find, METH_VARARGS, "Find test."},
    {"find_partial",  spats_find_partial, METH_VARARGS, "Find test."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
inittxspats(void)
{
    (void) Py_InitModule("txspats", SpatsMethods);
}
