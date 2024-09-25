#pragma once

// SPDX-FileCopyrightText: © 2020 Georg Sauthoff <mail@perfect_hash.tf>
// SPDX-License-Identifier: BSL-1.0

#include <stddef.h>
#include <cstdint>
#include <vector>


#include <exception>

namespace perfect_hash 
{

// popular general hash function
// originates from the sdbm package
// also used in GNU awk
// distributes better than djb2
// distributes worse than e.g. SipHash
// but is faster and good enough when the hash table size is a prime
static inline uint64_t hash_sdbm_64(const void *sP, size_t n, uint64_t param)
{
    uint64_t hash = 0;
    const unsigned char *s   = (const unsigned char*) sP;
    const unsigned char *end = s + n;

    // NB: 65599 is a prime that works well in practice
    // the parametrization isn't part of the original design,
    // i.e. this hash function wasn't designed as a parametized/seedable
    // hash function
    // however, this works good-enough to resolve collisions
    uint64_t k = 65599 + param;

    // NB: other versions of this function use some shifts
    //     however, modern optimizing compilers generate the same code,
    //     i.e. they replace the shifts with a multiplcation
    //     cf. https://godbolt.org/z/3Txh1n
    for (; s != end; ++s)
        hash = hash * k + *s;

    return hash;
}

// same as above but creates 32 bit hash values
// NB: on x86_64, saves some bytes on code size (e.g. 41 vs. 60 bytes with gcc 10.2)
//
// popular general hash function
// originates from the sdbm package
// also used in GNU awk
// distributes better than djb2
// distributes worse than e.g. SipHash
// but is faster and good enough when the hash table size is a prime
static inline uint32_t hash_sdbm_32(const void *sP, size_t n, uint32_t param)
{
    uint32_t hash = 0;
    const unsigned char *s   = (const unsigned char*) sP;
    const unsigned char *end = s + n;

    // NB: 65599 is a prime that works well in practice
    // the parametrization isn't part of the original design,
    // i.e. this hash function wasn't designed as a parametized/seedable
    // hash function
    // however, this works good-enough to resolve collisions
    uint32_t k = 65599 + param;

    // NB: other versions of this function use some shifts
    //     however, modern optimizing compilers generate the same code,
    //     i.e. they replace the shifts with a multiplcation
    //     cf. https://godbolt.org/z/3Txh1n
    for (; s != end; ++s)
        hash = hash * k + *s;

    return hash;
}


template< typename KEY, typename VALUE >
class table {
public:
    using data_entry_t = std::pair< KEY, VALUE >;
    using data_table_t = std::vector< data_entry_t >;

    using hash_function_t = uint32_t (*)(const void *p, size_t i, uint32_t param);

private:
    struct Gms_Phash_Bucket {
        uint32_t off;      // offset into table::idx_table
        uint8_t  n;        // secondary hash table size
        uint8_t  param;    // hash function parameter/seed to resolve collisions
    };

    std::vector< Gms_Phash_Bucket >     bkt_table;     // offsets into idx_table
    std::vector< uint32_t >             idx_table;
    uint32_t         bkt_table_n;
    uint32_t         idx_table_n;

    hash_function_t  hash_function_;
	uint32_t hash_for_build(const void *p, uint32_t i, uint32_t param) const
	{
		auto vector_data = static_cast< const data_entry_t * >( p );
		return hash_function_( &vector_data[ i ].first, sizeof(KEY), param);
	};
	uint32_t hash_for_lookup(const void *p, uint32_t i, uint32_t param) const 
	{
		(void)i;
		//const char *isin = p;
		return hash_function_(p, sizeof(KEY), param);
	};

public:
    using data_entry_t = std::pair< KEY, VALUE >;
    using data_table_t = std::vector< data_entry_t >;

    table() =default;
    table(const table &) =delete;
    table &operator=(const table &) =delete;
    table(table &&o) = default;
    table &operator=(table &&o) = default;

    table(const data_table_t& data, hash_function_t hash_function = hash_sdbm_32);
    inline uint32_t lookup(KEY key) const
    {
        uint32_t x = hash_for_lookup( &key, 0, 0);

        uint32_t i = (uint64_t)x * bkt_table_n >> 32;

        const Gms_Phash_Bucket *o = &bkt_table[ i ];

        // the expected case is that o->param is 0
        // if it's zero we could just assign x to y
        // thus, we could make this a conditional call
        // however, we don't do that in favor of equalizing the
        // runtime of the lookup for all elements (i.e. to minimize latency jitter)
        uint8_t y = hash_for_lookup( &key, 0, o->param);

        uint8_t j = (uint16_t)y * o->n >> 8;

        return idx_table[o->off + j];
    }
};

struct error : public std::exception
{
    int code {0};
    error(int code)
        : code(code)
    {
    }
    const char* what() const noexcept override {
        return "failed building perfect hash table";
    }
};

template< typename KEY, typename VALUE >
table< KEY, VALUE >::table( const data_table_t& data, hash_function_t hash_function )
    : hash_function_( hash_function )
{
    const void *p = data.data();
    std::size_t n = data.size();

    bkt_table_n = n/2;

    bkt_table.resize( bkt_table_n );

    std::vector< uint8_t > ns( bkt_table_n, 0 );

    for (uint32_t i = 0; i < n; ++i) {
        uint32_t x = hash_for_build(p, i, 0);
        uint32_t k = (uint64_t)x * bkt_table_n >> 32;
        ++ns[k];
        if (!ns[k]) {
            throw error(-2);
        }
    }

    std::vector< uint32_t* > vs( bkt_table_n, 0 );
    std::vector< uint32_t > ws( 2 * n, 0 );

    {
        auto t = ws.begin();
        for (uint32_t i = 0; i < bkt_table_n; ++i) {
            if (ns[i]) {
                vs[i] = &(*t);
                t += 2 * ns[i];
            }
        }
    }
    ns = std::vector< uint8_t >( bkt_table_n, 0 );

    for (uint32_t i = 0; i < n; ++i) {
        uint32_t x = hash_for_build(p, i, 0);
        uint32_t k = (uint64_t)x * bkt_table_n >> 32;
        vs[k][ns[k] * 2] = i;
        vs[k][ns[k] * 2 + 1] = x;
        ++ns[k];
    }

    bool col[256] = {0};
    uint32_t l = 0;
    for (uint32_t i = 0; i < bkt_table_n; ++i) {
        if (!ns[i])
            continue;
        if (ns[i] == 1) {
            bkt_table[i].off     = l;
            // already initialized to 0
            // NB: n = 0 or n = 1 is both fine
            // bkt_table[i].n     = 0;
            // bkt_table[i].param = 0;
            ++l;
        } else if (ns[i] > 1) {
            // printf("Collisions: %d\n", (int)ns[i]);
            // for (uint8_t k = 0; k < ns[i]; ++k) {
            //     printf("    %u -> %u\n", vs[i][k * 2], vs[i][k * 2 + 1]);
            // }

            uint32_t j = ns[i];
            bool done = true;
            uint8_t e = 0;
            for ( ; j < 256; ++j) {
                for (e = 0; e < 24; ++e) {
                    memset(col, 0, sizeof col);
                    done = true;
                    // printf("  Trying size %d (with param %d)\n", (int)j, (int)e);
                    for (uint8_t k = 0; k < ns[i]; ++k) {
                        uint8_t x = vs[i][k * 2 + 1];
                        if (e)
                            x = hash_for_build(p, vs[i][k * 2], e);
                        uint8_t  a = (uint16_t)x * j >> 8;
                        if (col[a]) {
                            done = false;
                            break;
                        }
                        col[a] = true;
                    }
                    if (done)
                        break;
                }
                if (done)
                    break;
            }
            if (!done) {
                throw error(-3);
            }
            bkt_table[i].off    = l;
            bkt_table[i].n      = j;
            bkt_table[i].param  = e;
            l += j;
        }
    }

    idx_table.resize( l );
 
    idx_table_n = l;
    for (uint32_t i = 0; i < bkt_table_n; ++i) {
        if (!ns[i])
            continue;
        for (uint8_t k = 0; k < ns[i]; ++k) {
            uint32_t a = vs[i][k * 2];
            uint32_t x = vs[i][k * 2 + 1];

            uint32_t j = (uint64_t)x * bkt_table_n >> 32;
            const Gms_Phash_Bucket *o = &bkt_table[ j ];

            uint8_t y = x;
            if (o->param)
                y = hash_for_build(p, vs[i][k * 2], o->param);
            
            uint8_t c = (uint16_t)y * o->n >> 8;
            idx_table[o->off + c] = a;
        }
    }
}

}

