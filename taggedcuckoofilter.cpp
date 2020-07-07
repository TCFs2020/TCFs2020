#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <fstream>
#include <math.h>
#include <bitset>
#include <x86intrin.h>

#include "taggedcuckoofilter.h"

using namespace std;

uint32_t FILTER_BIT_SIZE = (uint32_t) ((1LL << 20) * 32);
double FALSE_POSITIVE_RATE = 1E-6;
uint32_t BITS_PER_FINGERPRINT = 32;
uint32_t FILTER_MAX_ITEMS = (uint32_t) floor(CUCKOO_LOAD_FACTOR * FILTER_BIT_SIZE / BITS_PER_FINGERPRINT);
uint32_t FILTER_RAND_SEED = 1;

#define BITS_PER_UINT       32
#define BITS_PER_CHAR       8
static const uint8_t bit_mask[BITS_PER_CHAR] =
{
    0x01,  //00000001
    0x02,  //00000010
    0x04,  //00000100
    0x08,  //00001000
    0x10,  //00010000
    0x20,  //00100000
    0x40,  //01000000
    0x80   //10000000
};

/**
 *
 * This is the implementation of class TaggedCuckooFilter
 *
 */
TaggedCuckooFilter::TaggedCuckooFilter(uint32_t max_items)
{
    if(0 == max_items)
    {
        cout << "Error: Zero item in the set." << endl;
        exit(1);
    }

    num_items = 0;
    num_buckets = (uint32_t) floor(1.0 * FILTER_BIT_SIZE / (SLOTS_PER_BUCKET * BITS_PER_FINGERPRINT));
    num_buckets = (uint32_t) pow(2, ceil(log2(num_buckets)));
    fingerprint_bits = BITS_PER_FINGERPRINT;

    num_slots = num_buckets * SLOTS_PER_BUCKET;
    slot_array = new uint32_t [num_slots];
    memset(slot_array, 0, sizeof(uint32_t) * num_slots);
    filter_size = sizeof(uint32_t) * num_slots;

    item_hash_seed = 0;
    fingerprint_hash_seed = 0;
    GenerateItemFingerprintHashSeed();

    ResetMemoryAccesses();
}

TaggedCuckooFilter::~TaggedCuckooFilter()
{
    delete [] slot_array;
}

// generate an item hash seed and a fingerprint hash seed for the filter.
void TaggedCuckooFilter::GenerateItemFingerprintHashSeed()
{
    uint32_t salt_array[128] =
    {
        0xAAAAAAAA, 0x55555555, 0x33333333, 0xCCCCCCCC,
        0x66666666, 0x99999999, 0xB5B5B5B5, 0x4B4B4B4B,
        0xAA55AA55, 0x55335533, 0x33CC33CC, 0xCC66CC66,
        0x66996699, 0x99B599B5, 0xB54BB54B, 0x4BAA4BAA,
        0xAA33AA33, 0x55CC55CC, 0x33663366, 0xCC99CC99,
        0x66B566B5, 0x994B994B, 0xB5AAB5AA, 0xAAAAAA33,
        0x555555CC, 0x33333366, 0xCCCCCC99, 0x666666B5,
        0x9999994B, 0xB5B5B5AA, 0xFFFFFFFF, 0xFFFF0000,
        0xB823D5EB, 0xC1191CDF, 0xF623AEB3, 0xDB58499F,
        0xC8D42E70, 0xB173F616, 0xA91A5967, 0xDA427D63,
        0xB1E8A2EA, 0xF6C0D155, 0x4909FEA3, 0xA68CC6A7,
        0xC395E782, 0xA26057EB, 0x0CD5DA28, 0x467C5492,
        0xF15E6982, 0x61C6FAD3, 0x9615E352, 0x6E9E355A,
        0x689B563E, 0x0C9831A8, 0x6753C18B, 0xA622689B,
        0x8CA63C47, 0x42CC2884, 0x8E89919B, 0x6EDBD7D3,
        0x15B6796C, 0x1D6FDFE4, 0x63FF9092, 0xE7401432,
        0xEFFE9412, 0xAEAEDF79, 0x9F245A31, 0x83C136FC,
        0xC3DA4A8C, 0xA5112C8C, 0x5271F491, 0x9A948DAB,
        0xCEE59A8D, 0xB5F525AB, 0x59D13217, 0x24E7C331,
        0x697C2103, 0x84B0A460, 0x86156DA9, 0xAEF2AC68,
        0x23243DA5, 0x3F649643, 0x5FA495A8, 0x67710DF8,
        0x9A6C499E, 0xDCFB0227, 0x46A43433, 0x1832B07A,
        0xC46AFF3C, 0xB9C8FFF0, 0xC9500467, 0x34431BDF,
        0xB652432B, 0xE367F12B, 0x427F4C1B, 0x224C006E,
        0x2E7E5A89, 0x96F99AA5, 0x0BEB452A, 0x2FD87C39,
        0x74B2E1FB, 0x222EFD24, 0xF357F60C, 0x440FCB1E,
        0x8BBE030F, 0x6704DC29, 0x1144D12F, 0x948B1355,
        0x6D8FD7E9, 0x1C11A014, 0xADD1592F, 0xFB3C712E,
        0xFC77642F, 0xF9C4CE8C, 0x31312FB9, 0x08B0DD79,
        0x318FA6E7, 0xC040D23D, 0xC0589AA7, 0x0CA5C075,
        0xF874B172, 0x0CF914D5, 0x784D3280, 0x4E8CFEBC,
        0xC569F575, 0xCDB2A091, 0x2CC016B4, 0x5C5F4421
    };

    // select a rand item seed for hash functions.
    srand(FILTER_RAND_SEED);
    item_hash_seed = salt_array[(rand() % 128)];

    // select a rand fingerprint seed for hash functions.
    fingerprint_hash_seed = salt_array[(rand() % 128)];
}

// compute a hash beacon bucket index and a fingerprint key of an item.
void TaggedCuckooFilter::ComputeHashBeaconBucketFingerprintKey(uint64_t item_key, uint32_t & beacon_bucket_index, uint32_t & fingerprint_key)
{
    uint64_t hash_value = MurmurHash2_64((const char *)&item_key, sizeof(item_key), item_hash_seed);
    beacon_bucket_index = (uint32_t) ((hash_value >> 32) % num_buckets);
    fingerprint_key = (uint32_t) (hash_value % (1LL << fingerprint_bits));
    if(0 == fingerprint_key)
    {
        fingerprint_key = 1;
    }
}

// compute a hash alternative bucket index by XOR.
uint32_t TaggedCuckooFilter::ComputeHashAlternativeBucketIndex(uint32_t bucket_index, uint32_t fingerprint_key)
{
    uint32_t hash_value = MurmurHash3_32((const char *)&fingerprint_key, sizeof(fingerprint_key), fingerprint_hash_seed) % num_buckets;
    uint32_t alternative_bucket_index = (bucket_index ^ hash_value) % num_buckets;
    return alternative_bucket_index;
}

// store a fingerprint key into a bucket.
void TaggedCuckooFilter::StoreFingerprintKey(uint32_t bucket_index, uint32_t slot_index, uint32_t fingerprint_key)
{
    slot_array[(bucket_index * SLOTS_PER_BUCKET) + slot_index] = fingerprint_key;
}

// fetch a fingerprint key from a bucket.
void TaggedCuckooFilter::FetchFingerprintKey(uint32_t bucket_index, uint32_t slot_index, uint32_t & fingerprint_key)
{
    fingerprint_key = slot_array[(bucket_index * SLOTS_PER_BUCKET) + slot_index];
}

// check whether a bucket has a vacant slot.
bool TaggedCuckooFilter::IsBucketVacant(uint32_t bucket_index, uint32_t & slot_index)
{
    bool vacant_result = false;
    insert_accesses += 1;

    // check each slot in the bucket.
    slot_index = 0;
    for(size_t i = 0; i < SLOTS_PER_BUCKET; i++)
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(bucket_index, i, fingerprint_value);
        if(0 == fingerprint_value)
        {
            slot_index = i;
            vacant_result = true;
            break;
        }
    }

    return vacant_result;
}

// search a fingerprint key in the bucket.
bool TaggedCuckooFilter::SearchFingerprintKey(uint32_t bucket_index, uint32_t fingerprint_key)
{
    bool search_result = false;
    lookup_accesses += 1;

    // check each slot in the bucket.
    for(size_t i = 0; i < SLOTS_PER_BUCKET; i++)
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(bucket_index, i, fingerprint_value);
        if((0 != fingerprint_value) && (fingerprint_key == fingerprint_value))
        {
            search_result = true;
            break;
        }
    }

    return search_result;
}

// remove a fingerprint key from a bucket.
bool TaggedCuckooFilter::RemoveFingerprintKey(uint32_t bucket_index, uint32_t fingerprint_key)
{
    bool remove_result = false;
    delete_accesses += 1;

    // check each slot in the bucket.
    for(size_t i = 0; i < SLOTS_PER_BUCKET; i++)
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(bucket_index, i, fingerprint_value);
        if((0 != fingerprint_value) && (fingerprint_key == fingerprint_value))
        {
            StoreFingerprintKey(bucket_index, i, 0);
            remove_result = true;
            break;
        }
    }

    return remove_result;
}

// write the log of a bucket.
void TaggedCuckooFilter::WriteBucketLog(uint32_t bucket_index)
{
    ofstream bloom_result_txt;
    bloom_result_txt.open("results/bloom_result.txt", ios_base::app);
    bloom_result_txt << "Bucket [" << bucket_index << "]: " << endl;

    for(size_t i = 0; i < SLOTS_PER_BUCKET; i++)
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(bucket_index, i, fingerprint_value);
        bitset <32> char_bitset(fingerprint_value);
        bloom_result_txt << "Slot [" << i << "]: " << fingerprint_value << "(" << char_bitset << ")"<< endl;
    }

    bloom_result_txt.close();
}

// move items along the cuckoo path for inserting an item.
bool TaggedCuckooFilter::CuckooMove(uint32_t bucket_index, uint32_t fingerprint_key, uint32_t & num_kicks)
{
    bool move_result = false;

    uint32_t insert_bucket_index = bucket_index;
    uint32_t insert_fingerprint_key = fingerprint_key;
    uint32_t victim_bucket_index = 0;
    uint32_t victim_slot_index = 0;
    uint32_t victim_fingerprint_key = 0;

    while((false == move_result) && (num_kicks <= CUCKOO_MAX_KICKS))
    {
        // compute the alternative bucket index of the item.
        victim_bucket_index = ComputeHashAlternativeBucketIndex(insert_bucket_index, insert_fingerprint_key);

        if(true == IsBucketVacant(victim_bucket_index, victim_slot_index))
        {
            // store the insert fingerprint key into the victim bucket.
            StoreFingerprintKey(victim_bucket_index, victim_slot_index, insert_fingerprint_key);
            num_kicks++;
            move_result = true;
            break;
        }

        // kick out a random fingerprint key from the victim bucket.
        victim_slot_index = (uint32_t)(rand() % SLOTS_PER_BUCKET);
        FetchFingerprintKey(victim_bucket_index, victim_slot_index, victim_fingerprint_key);

        // store the fingerprint key into the victim bucket.
        StoreFingerprintKey(victim_bucket_index, victim_slot_index, insert_fingerprint_key);
        num_kicks++;

        insert_bucket_index = victim_bucket_index;
        insert_fingerprint_key = victim_fingerprint_key;
    }

    return move_result;
}

// insert an item into the filter.
bool TaggedCuckooFilter::Insert(uint64_t item_key)
{
    bool insert_result = false;

    // compute the beacon bucket index and the fingerprint key of the item.
    uint32_t beacon_bucket_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeHashBeaconBucketFingerprintKey(item_key, beacon_bucket_index, fingerprint_key);

    // check whether the beacon bucket has a vacant slot.
    uint32_t slot_index = 0;
    if(true == IsBucketVacant(beacon_bucket_index, slot_index))
    {
        // store the fingerprint key into the beacon bucket.
        StoreFingerprintKey(beacon_bucket_index, slot_index, fingerprint_key);
        num_items++;
        return true;
    }

    // compute the alternative bucket index of the item.
    uint32_t alternative_bucket_index = ComputeHashAlternativeBucketIndex(beacon_bucket_index, fingerprint_key);
    if(true == IsBucketVacant(alternative_bucket_index, slot_index))
    {
        // store the fingerprint key into the alternative bucket.
        StoreFingerprintKey(alternative_bucket_index, slot_index, fingerprint_key);
        num_items++;
        return true;
    }

    uint32_t * bucket_index_array = new uint32_t [CUCKOO_NUM_HASHES];
    bucket_index_array[0] = beacon_bucket_index;
    bucket_index_array[1] = alternative_bucket_index;

    // kick out other items to insert the fingerprint key.
    if(false == insert_result)
    {
        uint32_t num_kicks = 0;
        uint32_t rand_index = (uint32_t)(rand() % CUCKOO_NUM_HASHES);
        uint32_t victim_bucket_index = bucket_index_array[rand_index];
        uint32_t victim_slot_index = (uint32_t)(rand() % SLOTS_PER_BUCKET);
        uint32_t victim_fingerprint_key = 0;
        FetchFingerprintKey(victim_bucket_index, victim_slot_index, victim_fingerprint_key);

        // store the fingerprint key in the victim bucket.
        StoreFingerprintKey(victim_bucket_index, victim_slot_index, fingerprint_key);
        num_kicks++;

        // kick out other fingerprint keys.
        insert_result = CuckooMove(victim_bucket_index, victim_fingerprint_key, num_kicks);
        if(true == insert_result)
        {
            num_items++;
        }
    }

    delete [] bucket_index_array;
    return insert_result;
}

// lookup whether an item is in the filter.
bool TaggedCuckooFilter::Lookup(uint64_t item_key)
{
    bool lookup_result = false;

    // compute the beacon bucket index and the fingerprint key of the item.
    uint32_t beacon_bucket_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeHashBeaconBucketFingerprintKey(item_key, beacon_bucket_index, fingerprint_key);

    // search the fingerprint key in the beacon bucket.
    if(true == SearchFingerprintKey(beacon_bucket_index, fingerprint_key))
    {
        return true;
    }

    // search the fingerprint key in the alternative bucket.
    uint32_t alternative_bucket_index = ComputeHashAlternativeBucketIndex(beacon_bucket_index, fingerprint_key);
    if(true == SearchFingerprintKey(alternative_bucket_index, fingerprint_key))
    {
        return true;
    }

    return lookup_result;
}

// delete an item from the filter.
bool TaggedCuckooFilter::Delete(uint64_t item_key)
{
    bool delete_result = false;

    // compute the beacon bucket index and the fingerprint key of the item.
    uint32_t beacon_bucket_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeHashBeaconBucketFingerprintKey(item_key, beacon_bucket_index, fingerprint_key);

    // delete the fingerprint key from the beacon bucket.
    if(true == RemoveFingerprintKey(beacon_bucket_index, fingerprint_key))
    {
        num_items--;
        return true;
    }

    // delete the fingerprint key from the alternative bucket.
    uint32_t alternative_bucket_index = ComputeHashAlternativeBucketIndex(beacon_bucket_index, fingerprint_key);
    if(true == RemoveFingerprintKey(alternative_bucket_index, fingerprint_key))
    {
        num_items--;
        return true;
    }

    return delete_result;
}

// write the log of a cuckoo move path.
void TaggedCuckooFilter::WriteCuckooMoveLog(uint32_t bucket_index, uint64_t slot_index, uint32_t fingerprint_key, uint32_t num_kicks)
{
    ofstream cuckoo_log_txt;
    cuckoo_log_txt.open("results/bloom_result.txt", ios_base::app);

    bitset <32> char_bitset(fingerprint_key);
    cuckoo_log_txt << "Move " << num_kicks << ": [" << bucket_index << ", " << slot_index << ", "
                    << fingerprint_key << "(" << char_bitset << ")]" << endl;
    cuckoo_log_txt.close();
}

// calculate the theoretical size of the filter.
uint32_t TaggedCuckooFilter::CalculateFilterSize()
{
    return filter_size;
}

// calculate the bits per item of the filter.
double TaggedCuckooFilter::CalculateBitsPerItem()
{
    return (1.0 * filter_size * BITS_PER_CHAR / num_items);
}

// reset the operation memory accesses.
void TaggedCuckooFilter::ResetMemoryAccesses()
{
    insert_accesses = 0;
    lookup_accesses = 0;
    delete_accesses = 0;
}

// calculate the insert memory accesses.
uint64_t TaggedCuckooFilter::CalculateInsertAccesses()
{
    return insert_accesses;
}

// calculate the lookup memory accesses.
uint64_t TaggedCuckooFilter::CalculateLookupAccesses()
{
    return lookup_accesses;
}

// calculate the delete memory accesses.
uint64_t TaggedCuckooFilter::CalculateDeleteAccesses()
{
    return delete_accesses;
}

// write the log of the filter.
void TaggedCuckooFilter::WriteFilterLog()
{
    ofstream bloom_result_txt;
    bloom_result_txt.open("results/bloom_result.txt", ios_base::app);

    bloom_result_txt << "false_positive_rate: " << FALSE_POSITIVE_RATE
                     << ", load_factor: " << CUCKOO_LOAD_FACTOR
                     << ", num_items: " << num_items
                     << ", num_hashes: " << CUCKOO_NUM_HASHES
                     << ", slots_per_bucket: " << SLOTS_PER_BUCKET
                     << ", num_buckets: " << num_buckets
                     << ", fingerprint_bits: " << fingerprint_bits
                     << ", filter_size (bytes): " << CalculateFilterSize()
                     << ", bits_per_item: " << CalculateBitsPerItem() << endl;

    bloom_result_txt.close();
}
