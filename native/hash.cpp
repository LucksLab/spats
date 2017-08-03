/*
 * idea:
 * - given a set of N 64-bit inputs
 * - find a minimal-time perfect (or near) hash
 *
 * first pass:
 * (a) find any bits which don't contribute
 * (b) pack into low bits
 * (c) if minimal, then done; otherwise:
 * (d) for i in len(bits):
 *     - try A ^ (A >> i) -- no new collisions?
 *     - should be circular shift
 *     - once one is found, go back to step (b)
 */

#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


// reduce for debugging
#define BIT_LENGTH 10

typedef std::vector<uint64_t> Intlist;
typedef std::set<uint64_t> Intset;

inline uint64_t circular_shift(uint64_t n, int shift)
{
    if (shift >= 0) {
        uint64_t mask = ((1LL << shift) - 1);
        return ((n >> shift) | (n & mask) << (BIT_LENGTH - shift));
    }
    else {
        uint64_t mask = ~((1LL << (BIT_LENGTH + shift)) - 1);
        return ((n << -shift) | (n & mask) >> (BIT_LENGTH + shift));
    }
}

inline std::string sformat(const char * formatStr, ...) 
{
    /* The following is derived from one of the answers at:
     *  http://stackoverflow.com/questions/2342162/stdstring-formating-like-sprintf  */
    int size = 64;
    std::string result;
    va_list ap;
    while (true)
    {
        result.resize(size);
        va_start(ap, formatStr);
        int n = vsnprintf(const_cast<char *>(result.c_str()), size, formatStr, ap);
        va_end(ap);
        if (n > -1  &&  n < size)
            return result.substr(0, n);
        size = (n > -1 ?  n + 1  :  size * 2);
    } // while
}


std::string
bitstring(uint64_t n)
{
    char bitbuf[BIT_LENGTH + 1];
    uint64_t mask = 1LL;
    for (int i = 0; i < BIT_LENGTH; ++i) {
        bitbuf[BIT_LENGTH - 1 - i] = ((n & mask) ? '1' : '0');
        mask = (mask << 1);
    }
    bitbuf[BIT_LENGTH] = 0;
    return std::string(bitbuf);
}
#define BSTR(n) (bitstring(n).c_str())

#define NUM_STATE_REGS 4
class State
{
public:
    uint64_t r[NUM_STATE_REGS];
    State() { reset(); }
    void reset() { for (int i = 0; i < NUM_STATE_REGS; ++i) r[i] = 0LL; }
    std::string dump() { return sformat("0/1: %s %s\n2/3:%s %s", BSTR(r[0]), BSTR(r[1]), BSTR(r[2]), BSTR(r[3])); }
};

class Op
{
public:
    virtual ~Op() { }
    virtual void perform(State * state) = 0;
    virtual std::string make_asm() { return "OPNYI"; }
    virtual std::string desc() { return "DESCNYI"; }
};

class Mask :
    public virtual Op
{
public:
    int src_reg;
    int dst_reg;
    uint64_t mask;
    Mask(int src, int dst, uint64_t m) : src_reg(src), dst_reg(dst), mask(m) { }
    virtual void perform(State * state) { state->r[dst_reg] = (state->r[src_reg] & mask); }
    virtual std::string desc() { return sformat("r[%d] = (r[%d] & %s)", dst_reg, src_reg, BSTR(mask)); }
    virtual ~Mask() { }
};

class Or :
    public virtual Op
{
public:
    int src1_reg;
    int src2_reg;
    int dst_reg;
    Or(int src1, int src2, int dst) : src1_reg(src1), src2_reg(src2), dst_reg(dst) { }
    virtual void perform(State * state) { state->r[dst_reg] = (state->r[src1_reg] | state->r[src2_reg]); }
    virtual std::string desc() { return sformat("r[%d] = (r[%d] | r[%d])", dst_reg, src1_reg, src2_reg); }
    virtual ~Or() { }
};

class Xor :
    public virtual Op
{
public:
    int src1_reg;
    int src2_reg;
    int dst_reg;
    Xor(int src1, int src2, int dst) : src1_reg(src1), src2_reg(src2), dst_reg(dst) { }
    virtual void perform(State * state) { state->r[dst_reg] = (state->r[src1_reg] ^ state->r[src2_reg]); }
    virtual std::string desc() { return sformat("r[%d] = (r[%d] ^ r[%d])", dst_reg, src1_reg, src2_reg); }
    virtual ~Xor() { }
};

class Shift :
    public virtual Op
{
public:
    int src_reg;
    int dst_reg;
    int shift;
    Shift(int src, int dst, int s) : src_reg(src), dst_reg(dst), shift(s) { }
    virtual void perform(State * state) { state->r[dst_reg] = (shift >= 0) ? (state->r[src_reg] >> shift) : (state->r[src_reg] << -shift); }
    virtual std::string desc() { return sformat("r[%d] = (r[%d] %s %d)", dst_reg, src_reg, (shift >= 0 ? ">>" : "<<"), (shift >= 0 ? shift : -shift)); }
    virtual ~Shift() { }
};

class CShift :
    public virtual Op
{
public:
    int src_reg;
    int dst_reg;
    int shift;
    CShift(int src, int dst, int s) : src_reg(src), dst_reg(dst), shift(s) { }
    virtual void perform(State * state) { state->r[dst_reg] = circular_shift(state->r[src_reg], shift); }
    virtual std::string desc() { return sformat("r[%d] = (r[%d] %s %d)", dst_reg, src_reg, (shift >= 0 ? ">>]" : "[<<"), (shift >= 0 ? shift : -shift)); }
    virtual ~CShift() { }
};

class Fn
{
private:
    std::vector<Op *> _ops;
public:
    void addOp(Op * op) { _ops.push_back(op); }
    void perform(State * state) { for (std::vector<Op *>::iterator it = _ops.begin(); it != _ops.end(); ++it) (*it)->perform(state); }
    std::string desc()
    {
        std::string desc(sformat("F[%d]:\n", _ops.size()));
        for (std::vector<Op *>::iterator it = _ops.begin(); it != _ops.end(); ++it)
            desc += sformat("  %s\n", (*it)->desc().c_str());
        return desc + "------\n";
    }
};

void
find_used_bits(Intlist * inputs, Intlist * used_bits)
{
    uint64_t mask = 0LL;

    for (uint64_t i = 0; i < BIT_LENGTH; ++i) {
        uint64_t cur_mask = ((1LL << i) | mask);
        //printf("    %s:\n", BSTR(cur_mask));
        Intset s;
        for (Intlist::iterator it = inputs->begin(); it != inputs->end(); ++it) {
            //printf("       %s\n", BSTR((*it)&(~cur_mask)));
            s.insert((*it) & (~cur_mask));
        }
        if (s.size() < inputs->size()) {
            //printf("    (used)\n");
            used_bits->push_back(i);
        }
        else {
            mask = cur_mask;
            //printf("    **unused\n");
        }
    }
}

class Range
{
public:
    uint64_t loc;
    uint64_t len;
    Range() : loc(0LL), len(0LL) { }
    Range(uint64_t _loc, uint64_t _len) : loc(_loc), len(_len) { }
    uint64_t end() { return loc + len; }
    void reset() { loc = len = 0LL; }
};

void
find_holes(Intlist * used_bits, std::vector<Range> * holes)
{
    Range r;
    for (uint64_t i = 0LL; i < BIT_LENGTH; ++i) {
        if (used_bits->end() != std::find(used_bits->begin(), used_bits->end(), i)) {
            if (r.len) {
                holes->push_back(r);
                r.reset();
            }
            continue;
        }
        if (!r.len) {
            r.loc = i;
            r.len = 1;
        }
        else if (r.end() == i) {
            ++r.len;
        }
        else {
            assert(0);
        }
    }
    if (0 < r.len)
        holes->push_back(r);
}

int
pack_bits(Intlist * inputs, Fn * fn)
{
    Intlist used_bits;
    find_used_bits(inputs, &used_bits);
    std::vector<Range> holes;
    find_holes(&used_bits, &holes);
    int num_shifted_so_far = 0;
    for (std::vector<Range>::iterator it = holes.begin(); it != holes.end(); ++it) {
        Range hole(*it);
        printf(" Hole: %d:%d\n", (int)hole.loc, (int)hole.len);
        /* the top hole doesn't require any work */
        if (hole.end() == BIT_LENGTH)
            continue;
        hole.loc -= num_shifted_so_far;
        if (0 == hole.loc) {
            /* special case of simple right shift */
            fn->addOp(new Shift(0, 0, (int)hole.len));
        }
        else {
            /* op: ((val & highmask) >> hole.len) | (val & lowmask) */
            uint64_t lmask = ((1LL << hole.loc) - 1LL);
            uint64_t hmask = ~((1LL << hole.end()) - 1LL);
            //printf("  lm: %s, hm: %s\n", BSTR(lmask), BSTR(hmask));
            fn->addOp(new Mask(0, 1, hmask));
            fn->addOp(new Mask(0, 2, lmask));
            fn->addOp(new Shift(1, 1, (int)hole.len));
            fn->addOp(new Or(1, 2, 0));
        }
        num_shifted_so_far += hole.len;
    }
    return (BIT_LENGTH - num_shifted_so_far);
}

void
show_inputs(Intlist * ins, bool binary = false)
{
    for (Intlist::iterator it = ins->begin(); it != ins->end(); ++it) {
        if (binary)
            printf("%s  %d\n", BSTR(*it), (int)*it);
        else
            printf("%d ", (int)*it);
    }
    if (!binary)
        printf("\n");
}

void
map_inputs(Intlist * ins, Intlist * mapped, Fn * f)
{
    State s;
    for (Intlist::iterator it = ins->begin(); it != ins->end(); ++it) {
        s.reset();
        s.r[0] = *it;
        f->perform(&s);
        mapped->push_back(s.r[0]);
    }
}

uint64_t
max_int(Intlist * ins)
{
    uint64_t max = 0;
    for (Intlist::iterator it = ins->begin(); it != ins->end(); ++it)
        if (*it > max)
            max = *it;
    return max;
}

void
find_xor(Intlist * inputs, Fn *f, int k)
{
    // TAI...
    f->addOp(new CShift(0, 1, (1 + (k % 4)) * (k & 1 ? 1 : -1)));
    f->addOp(new Xor(0, 1, 0));
}

bool
find_fn(Intlist * inputs, Fn *f)
{
    Intlist cur(*inputs);

    int tries = 0;

    while (true)
    {
        printf("Iteration %d...\n", tries);
        int bits = pack_bits(&cur, f);
        printf("  bits: %d\n", bits);
        Intlist mapped;
        map_inputs(inputs, &mapped, f);
        if (max_int(&mapped) <= mapped.size()) {
            printf(" ** Optimal!\n");
            return true;
        }
        if (max_int(&mapped) <= mapped.size() * 2) {
            printf(" ** Close to optimal...\n");
            return true;
        }
        if (++tries > 20) {
            printf("  aborting...\n");
            return false;
        }
        find_xor(&cur, f, tries);
        cur.clear();
        map_inputs(inputs, &cur, f);
        //show_inputs(&cur, true);
    }
}

void
make_ins(uint64_t * nums, Intlist * ins)
{
    for(int i = 0; i == 0  ||  0 != nums[i]; ++i)
        ins->push_back(nums[i]);
}

uint64_t set_1[] = { 1LL, 3LL, 5LL, 9LL, 16LL, 17LL, 19LL, 21LL, 0LL };
uint64_t set_2[] = { 0LL, 1LL, 2LL, 4LL, 8LL, 16LL };
uint64_t * 
set_3(int k)
{
    //sranddev();
    uint64_t * res = (uint64_t *)malloc((k + 1) * sizeof(uint64_t));
    for (int i = 0; i < k; ++i)
        res[i] = (unsigned long long)((rand() % 1024));
    res[k] = 0LL;
    return res;
}

void
thash()
{
    Intlist ins;
    make_ins((uint64_t *)set_3(34), &ins);
    show_inputs(&ins, true);
    printf("------\n");

    Fn f;
    if(!find_fn(&ins, &f))
        return;
    printf("%s", f.desc().c_str());

    Intlist mapped;
    map_inputs(&ins, &mapped, &f);
    show_inputs(&mapped, true);
}
