#ifndef TAGGEDCUCKOOFILTER_H_INCLUDED
#define TAGGEDCUCKOOFILTER_H_INCLUDED

#include "murmurhash3.h"

extern uint32_t FILTER_BIT_SIZE;
extern double   FALSE_POSITIVE_RATE;
extern uint32_t BITS_PER_FINGERPRINT;
extern uint32_t FILTER_MAX_ITEMS;
extern uint32_t FILTER_RAND_SEED;

#define CUCKOO_LOAD_FACTOR      0.95  // 0.95
#define CUCKOO_NUM_HASHES       2  // k=2
#define SLOTS_PER_BUCKET        4  // b=4
#define CUCKOO_MAX_KICKS        500

class TaggedCuckooFilter
{
public:
    uint32_t num_items;
    uint32_t num_buckets;
    uint32_t item_hash_seed;
    uint32_t fingerprint_hash_seed;
    uint32_t fingerprint_bits;
    uint32_t num_slots;
    uint32_t * slot_array;
    uint32_t filter_size;
    uint64_t insert_accesses;
    uint64_t lookup_accesses;
    uint64_t delete_accesses;

private:
    void GenerateItemFingerprintHashSeed();
    void ComputeHashBeaconBucketFingerprintKey(uint64_t item_key, uint32_t & beacon_bucket_index, uint32_t & fingerprint_key);
    uint32_t ComputeHashAlternativeBucketIndex(uint32_t bucket_index, uint32_t fingerprint_key);
    void StoreFingerprintKey(uint32_t bucket_index, uint32_t slot_index, uint32_t fingerprint_key);
    void FetchFingerprintKey(uint32_t bucket_index, uint32_t slot_index, uint32_t & fingerprint_key);
    bool IsBucketVacant(uint32_t bucket_index, uint32_t & slot_index);
    bool SearchFingerprintKey(uint32_t bucket_index, uint32_t fingerprint_key);
    bool RemoveFingerprintKey(uint32_t bucket_index, uint32_t fingerprint_key);
    void WriteBucketLog(uint32_t bucket_index);
    bool CuckooMove(uint32_t bucket_index, uint32_t fingerprint_key, uint32_t & num_kicks);
    void WriteCuckooMoveLog(uint32_t bucket_index, uint64_t slot_index, uint32_t fingerprint_key, uint32_t num_kicks);

public:
    TaggedCuckooFilter(uint32_t max_items);
    ~TaggedCuckooFilter();
    bool Insert(uint64_t item_key);
    bool Lookup(uint64_t item_key);
    bool Delete(uint64_t item_key);
    uint32_t CalculateFilterSize();
    double CalculateBitsPerItem();
    void ResetMemoryAccesses();
    uint64_t CalculateInsertAccesses();
    uint64_t CalculateLookupAccesses();
    uint64_t CalculateDeleteAccesses();
    void WriteFilterLog();
};

#endif // TAGGEDCUCKOOFILTER_H_INCLUDED
