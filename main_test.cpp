#include <x86intrin.h>

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <algorithm>

#include "rans_byte_test.h"
#include "platform.h"

template<class T>
T saveDiv(T a, T b) {
    return (a + b - 1) / b;
}

// This is just the sample program. All the meat is in rans_byte.h.

static void panic(const char *fmt, ...)
{
    va_list arg;

    va_start(arg, fmt);
    fputs("Error: ", stderr);
    vfprintf(stderr, fmt, arg);
    va_end(arg);
    fputs("\n", stderr);

    exit(1);
}

static uint8_t* read_file(char const* filename, size_t* out_size)
{
    FILE* f = fopen(filename, "rb");
    if (!f)
        panic("file not found: %s\n", filename);

    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    uint8_t* buf = new uint8_t[size];
    if (fread(buf, size, 1, f) != 1)
        panic("read failed\n");

    fclose(f);
    if (out_size)
        *out_size = size;

    return buf;
}

// ---- Stats

struct SymbolStats
{
    uint32_t freqs[256];
    uint32_t cum_freqs[257];

    void count_freqs(uint8_t const* in, size_t nbytes);
    void calc_cum_freqs();
    void normalize_freqs(uint32_t target_total);
};

void SymbolStats::count_freqs(uint8_t const* in, size_t nbytes)
{
    for (int i=0; i < 256; i++)
        freqs[i] = 0;

    for (size_t i=0; i < nbytes; i++)
        freqs[in[i]]++;
}

void SymbolStats::calc_cum_freqs()
{
    cum_freqs[0] = 0;
    for (int i=0; i < 256; i++)
        cum_freqs[i+1] = cum_freqs[i] + freqs[i];
}

void SymbolStats::normalize_freqs(uint32_t target_total)
{
    assert(target_total >= 256);
    
    calc_cum_freqs();
    uint32_t cur_total = cum_freqs[256];
    
    // resample distribution based on cumulative freqs
    for (int i = 1; i <= 256; i++)
        cum_freqs[i] = ((uint64_t)target_total * cum_freqs[i])/cur_total;

    // if we nuked any non-0 frequency symbol to 0, we need to steal
    // the range to make the frequency nonzero from elsewhere.
    //
    // this is not at all optimal, i'm just doing the first thing that comes to mind.
    for (int i=0; i < 256; i++) {
        if (freqs[i] && cum_freqs[i+1] == cum_freqs[i]) {
            // symbol i was set to zero freq

            // find best symbol to steal frequency from (try to steal from low-freq ones)
            uint32_t best_freq = ~0u;
            int best_steal = -1;
            for (int j=0; j < 256; j++) {
                uint32_t freq = cum_freqs[j+1] - cum_freqs[j];
                if (freq > 1 && freq < best_freq) {
                    best_freq = freq;
                    best_steal = j;
                }
            }
            assert(best_steal != -1);

            // and steal from it!
            if (best_steal < i) {
                for (int j = best_steal + 1; j <= i; j++)
                    cum_freqs[j]--;
            } else {
                assert(best_steal > i);
                for (int j = i + 1; j <= best_steal; j++)
                    cum_freqs[j]++;
            }
        }
    }

    // calculate updated freqs and make sure we didn't screw anything up
    assert(cum_freqs[0] == 0 && cum_freqs[256] == target_total);
    for (int i=0; i < 256; i++) {
        if (freqs[i] == 0)
            assert(cum_freqs[i+1] == cum_freqs[i]);
        else
            assert(cum_freqs[i+1] > cum_freqs[i]);

        // calc updated freq
        freqs[i] = cum_freqs[i+1] - cum_freqs[i];
    }
}

int main()
{
    size_t in_size;
    uint8_t* in_bytes = read_file("book11", &in_size);

    static const uint32_t prob_bits = 16;
    static const uint32_t prob_scale = 1 << prob_bits;

    SymbolStats stats;
    stats.count_freqs(in_bytes, in_size);
    stats.normalize_freqs(prob_scale);

    // cumlative->symbol table
    // this is super brute force
    uint32_t cum2sym[prob_scale];
    for (int s=0; s < 256; s++)
        for (uint32_t i=stats.cum_freqs[s]; i < stats.cum_freqs[s+1]; i++)
            cum2sym[i] = s;

    static size_t out_max_size = 32<<22;
    uint32_t* out_buf = new uint32_t[out_max_size];
    uint32_t* state_buf = new uint32_t[out_max_size];
    uint8_t* dec_bytes = new uint8_t[in_size];

    // try rANS encode
    uint32_t *rans_begin, *state_begin;
    RansEncSymbol esyms[256];
    RansDecSymbol dsyms[256];

    for (int i=0; i < 256; i++) {
        RansEncSymbolInit(&esyms[i], stats.cum_freqs[i], stats.freqs[i], prob_bits);
        RansDecSymbolInit(&dsyms[i], stats.cum_freqs[i], stats.freqs[i]);
    }

    // ---- regular rANS encode/decode. Typical usage.

    printf("rANS encode:\n");
    for (int run=0; run < 1; run++) {
        RansState rans;
        RansEncInit(&rans);

        uint32_t* ptr = out_buf + out_max_size; // *end* of output buffer
        uint32_t* xptr = state_buf + out_max_size;
        for (size_t i = in_size; i > 0; i--) { // NB: working in reverse!
            int s = in_bytes[i - 1];
            //printf("Encode %c\n", s);
            RansEncPutSymbol(&rans, &ptr, &xptr, &esyms[s]);
        }
        RansEncFlush(&rans, &ptr);
        rans_begin = ptr;
        state_begin = xptr;
    }
    printf("rANS: %d bytes\n", (int) (out_buf + out_max_size - rans_begin));

    memset(dec_bytes, 0xcc, in_size);

    double start_time = timer();
    // try normal rANS decode
    for (int run=0; run < 1; run++) {
        RansState rans;
        uint32_t* ptr = rans_begin;
        RansDecInit(&rans, &ptr);

        for (size_t i=0; i < in_size; i++) {
            uint32_t s = cum2sym[RansDecGet(&rans, prob_bits)];
            dec_bytes[i] = (uint8_t) s;
            RansDecAdvanceSymbol(&rans, &ptr, &dsyms[s], prob_bits);
            //printf("Symbol is: %c\n", s);
        }
    }
    double dec_time = timer() - start_time;
    printf("dec time with preprocessed rans: %.5f s\n", dec_time);

    // check decode results
    if (memcmp(in_bytes, dec_bytes, in_size) == 0)
        printf("decode ok!\n");
    else
        printf("ERROR: bad decoder!\n");

    memset(dec_bytes, 0xcc, in_size);

    start_time = timer();
    // try non preprocessed normal rANS decode
    for (int run=0; run < 1; run++) {
        RansState rans;
        uint32_t* ptr = rans_begin;
        RansDecInit(&rans, &ptr);

        for (size_t i=0; i < in_size; i++) {
            uint32_t cum_freq = RansDecGet(&rans, prob_bits);
            uint32_t s = std::find_if(stats.cum_freqs, stats.cum_freqs + 257, [cum_freq](int v) { return v > cum_freq; })
                    - stats.cum_freqs - 1;

            dec_bytes[i] = (uint8_t) s;
            RansDecAdvance(&rans, &ptr, stats.cum_freqs[s], stats.freqs[s], prob_bits);
        }
    }
    dec_time = timer() - start_time;
    printf("dec time with non preprocessed rans: %.5f s\n", dec_time);

    // check decode results
    if (memcmp(in_bytes, dec_bytes, in_size) == 0)
        printf("decode ok!\n");
    else
        printf("ERROR: bad decoder!\n");

    memset(dec_bytes, 0xcc, in_size);

    const int FINAL_STATE_SIZE = 1;
    int bitstream_n = out_buf + out_max_size - rans_begin - FINAL_STATE_SIZE;
    auto split_size = saveDiv<size_t>(bitstream_n, 8);

    uint8_t* dec_bytes_split[8];
    for (int i = 0; i < 8; i++) {
        dec_bytes_split[i] = new uint8_t[in_size];
    }

    RansState rans[8] __attribute__((aligned(32)));
    uint32_t* ptr[8] __attribute__((aligned(32))) = { rans_begin };
    int32_t offset[8] __attribute__((aligned(32))) = { FINAL_STATE_SIZE };
    int32_t offsetBounds[8] __attribute__((aligned(32)));

    RansDecInit(&rans[0], &ptr[0]);
    for (int i = 1; i < 8; i++) {
        ptr[i] = rans_begin + FINAL_STATE_SIZE + split_size * i;
        rans[i] = *(state_begin + split_size * i);
        RansDecInitAdvance(&rans[i], &ptr[i]);

        offset[i] = ptr[i] - rans_begin;
        offsetBounds[i - 1] = offset[i] - 1;
    }

    offsetBounds[7] = out_buf + out_max_size - rans_begin;

    __m256i probMask = _mm256_set1_epi32(prob_scale - 1); // (1 << prob_bits) - 1
    __m256i offsetBoundsSimd = _mm256_load_si256(reinterpret_cast<__m256i*>(offsetBounds));

    __m256i ransSimd = _mm256_load_si256(reinterpret_cast<__m256i*>(rans));
    __m256i offsetSimd = _mm256_load_si256(reinterpret_cast<__m256i*>(offset));
    __m256i endIndex = _mm256_setzero_si256(); // 0
    __m256i runningFlag = _mm256_set1_epi64x(-1); // all bits 1

    start_time = timer();
    for (int symbolCount = 0; ; symbolCount++) {
        // Extract Symbol

        // currSymbolCdf = ransSimd & symbolMask;
        __m256i currSymbolCdf = _mm256_and_si256(ransSimd, probMask);

        // currSymbol = cum2sym[currSymbolCdf];
        __m256i currSymbol = _mm256_i32gather_epi32(reinterpret_cast<const int*>(cum2sym), currSymbolCdf, 4);

        /*
        __m256i currSymbol = _mm256_setzero_si256();
        for (int symbol = 1; symbol < sizeof(stats.cum_freqs); symbol++) {
            // currCdfValue = stats.cum_freqs[symbol];
            //__m256i currCdfValue = _mm256_i32gather_epi32(reinterpret_cast<const int*>(stats.cum_freqs), _mm256_set1_epi32(symbol), 4);
            __m256i currCdfValue = _mm256_set1_epi32(stats.cum_freqs[symbol]);

            // flags = currCdfValue > currSymbolCdf;
            __m256i flags = _mm256_cmpgt_epi32(currCdfValue, currSymbolCdf);

            // newlyFound = flags - foundFlags;
            __m256i newlyFound = _mm256_and_si256(flags, _mm256_cmpeq_epi32(currSymbol, _mm256_setzero_si256()));

            if (_mm256_movemask_epi8(newlyFound)) {
                // currSymbol = newlyFound ? symbol - 1 : currSymbol
                currSymbol = _mm256_blendv_epi8(currSymbol, _mm256_set1_epi32(symbol), newlyFound);

                if (!~_mm256_movemask_epi8(flags)) break;
            }
        }

        // currSymbol -= 1
        currSymbol = _mm256_sub_epi32(currSymbol, _mm256_set1_epi32(1));*/

        /*uint32_t cum_freqs[8] __attribute__((aligned(32)));
        _mm256_store_si256((__m256i*)cum_freqs, currSymbolCdf);
        uint32_t s[8] __attribute__((aligned(32)));
        for (int i = 0; i < 8; i++) {
            s[i] = std::find_if(stats.cum_freqs, stats.cum_freqs + 257, [&](int v) { return v > cum_freqs[i]; }) - stats.cum_freqs - 1;
        }
        currSymbol = _mm256_load_si256(reinterpret_cast<__m256i*>(s));*/

        // TODO: store output properly

        int out[8] __attribute__((aligned(256)));
        _mm256_store_si256((__m256i*)out, currSymbol);
        for (int j = 0; j < 8; j++) {
            dec_bytes_split[j][symbolCount] = out[j];
        }

        // Advance states

        // currSymbolFreq = stats.freqs[currSymbol];
        // currSymbolCleanCdf = stats.cum_freqs[currSymbol];
        __m256i currSymbolFreq = _mm256_i32gather_epi32(reinterpret_cast<const int*>(stats.freqs), currSymbol, 4);
        __m256i currSymbolCleanCdf = _mm256_i32gather_epi32(reinterpret_cast<const int*>(stats.cum_freqs), currSymbol, 4);

        // ransSimd = (ransSimd >> prob_bits) * currSymbolFreq + currSymbolCdf - currSymbolCleanCdf;
        ransSimd = _mm256_mullo_epi32(_mm256_srli_epi32(ransSimd, prob_bits), currSymbolFreq);
        ransSimd = _mm256_sub_epi32(_mm256_add_epi32(ransSimd, currSymbolCdf), currSymbolCleanCdf);

        // Check renormalization

        // renormMask = ransSimd < RANS_BYTE_L; /* hack for unsigned comparison */
        __m256i renormMask = _mm256_xor_si256(ransSimd, _mm256_set1_epi32(0x80000000));
        renormMask = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L - 0x80000000),renormMask);

        // renormMask &= runningFlag;
        renormMask = _mm256_and_si256(renormMask, runningFlag);

        if (_mm256_movemask_epi8(renormMask)) {
            // nextRenormStates = rans_begin[offsetSimd];
            __m256i nextRenormStates = _mm256_i32gather_epi32(reinterpret_cast<const int*>(rans_begin), offsetSimd, 4);

            // renormedRans = ransSimd << BITS_WRITEOUT | nextRenormStates;
            __m256i renormedRans = _mm256_or_si256(_mm256_slli_epi32(ransSimd, BITS_WRITEOUT), nextRenormStates);

            // ransSimd = renormMask ? renormedRans : ransSimd;
            ransSimd = _mm256_blendv_epi8(ransSimd, renormedRans, renormMask);

            // threadsEnded = (offsetSimd == offsetBoundsSimd) & renormMask
            __m256i threadsEnded = _mm256_and_si256(_mm256_cmpeq_epi32(offsetSimd, offsetBoundsSimd), renormMask);

            if (_mm256_movemask_epi8(threadsEnded)) {
                // endIndex = threadsEnded ? symbolCount : outLength;
                endIndex = _mm256_blendv_epi8(endIndex, _mm256_set1_epi32(symbolCount), threadsEnded);

                // runningFlag &= ~threadsEnded;
                runningFlag = _mm256_and_si256(runningFlag, ~threadsEnded);
            }

            // offsetSimd += renormMask & 0x01
            offsetSimd = _mm256_add_epi32(offsetSimd, _mm256_and_si256(renormMask, _mm256_set1_epi32(0x01)));

            if (!_mm256_movemask_epi8(runningFlag)) {
                break;
            }
        }
    }
    dec_time = timer() - start_time;
    printf("dec time with simd: %.5f s\n", dec_time);

    int endIndexes[8] __attribute__((aligned(256)));
    _mm256_store_si256((__m256i*)endIndexes, endIndex);

    uint8_t *dec_bytes_ptr = dec_bytes;
    for (int i = 0; i < 8; i++) {
        int len = endIndexes[i] + 1;
        memcpy(dec_bytes_ptr, dec_bytes_split[i], len);
        dec_bytes_ptr += len;
    }

    // check decode results
    if (memcmp(in_bytes, dec_bytes, in_size) == 0)
        printf("simd decode ok!\n");
    else
        printf("ERROR: bad simd decoder!\n");

    delete[] out_buf;
    delete[] dec_bytes;
    delete[] in_bytes;
    return 0;
}
