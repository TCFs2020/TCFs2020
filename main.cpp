#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <dir.h>
#include <math.h>
#include <string.h>

#include "taggedcuckoofilter.h"

using namespace std;

#define MAX_SIZE_RUNS       3  // 32Mb, 128Mb, 512Mb.
#define MAX_TEST_RUNS       10  // 10 test runs.
#define MAX_FPR_RUNS        3  // 1E-2, 1E-4, 1E-6.
#define SLOT_BIT_SIZE       32  // 32 bits per slot.
uint32_t SLOT_INDEX_BITS = 20;  // 2^20=1M, 2^22=4M, 2^24=16M.
uint64_t START_ITEM_KEY = 1;

// generate an item array.
void GenerateItemArray(uint64_t item_array[], uint64_t max_items)
{
    START_ITEM_KEY = (((uint64_t)rand()) << 32) + (uint64_t)rand();

    for(uint64_t i  = 0; i < max_items; i++)
    {
        item_array[i] = (uint64_t)(START_ITEM_KEY + i);
    }
}

// generate a test array.
void GenerateTestArray(uint64_t item_array[], uint64_t max_items, uint64_t test_array[], uint64_t max_tests, uint32_t true_test_fraction)
{
    uint64_t true_tests = (uint64_t) (1.0 * max_tests * true_test_fraction / 100);
    uint64_t false_tests = max_tests - true_tests;

    for(uint64_t i = 0; i < true_tests; i++)
    {
        test_array[i] = item_array[i % max_items];
    }

    for(uint64_t j = 0; j < false_tests; j++)
    {
        test_array[true_tests + j] = item_array[max_items - 1] + j + 1;
    }
}

// test the bloom filter.
void TestBloomFilter()
{
    // create the result files.
    mkdir("results");
    remove("results/bloom_result.txt");
    remove("results/insert_result.txt");
    remove("results/lookup_result.txt");
    remove("results/delete_result.txt");

    ofstream insert_result_txt;
    insert_result_txt.open("results/insert_result.txt", ios_base::app);
    ofstream lookup_result_txt;
    lookup_result_txt.open("results/lookup_result.txt", ios_base::app);
    ofstream delete_result_txt;
    delete_result_txt.open("results/delete_result.txt", ios_base::app);

    for(size_t h = 0; h < MAX_SIZE_RUNS; h++)  // 32Mb, 128Mb, 512Mb.
    {
        FALSE_POSITIVE_RATE = 1E-2;
        for(size_t i = 0; i < MAX_FPR_RUNS; i++)  // 1E-2, 1E-4, 1E-6.
        {
			FILTER_BIT_SIZE = (uint32_t) ((1LL << SLOT_INDEX_BITS) * SLOT_BIT_SIZE);
			FILTER_MAX_ITEMS = (uint32_t) floor(CUCKOO_LOAD_FACTOR * FILTER_BIT_SIZE / BITS_PER_FINGERPRINT);

            uint32_t table_occupancy_fraction = 10;
            while(table_occupancy_fraction <= 100)  // 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 100%.
            {
                uint64_t max_items = (uint64_t) (1.0 * FILTER_MAX_ITEMS * table_occupancy_fraction / 100);
                uint64_t max_tests = 10 * max_items;

                for(size_t j = 0; j < MAX_TEST_RUNS; j++)  // 10 test runs.
                {
                    // generate an item array.
                    uint64_t * item_array = new uint64_t [max_items];
                    memset(item_array, 0, sizeof(uint64_t) * max_items);
                    GenerateItemArray(item_array, max_items);

                    // create a bloom filter.
                    class TaggedCuckooFilter * bloom_filter = new class TaggedCuckooFilter(max_items);
                    FILTER_RAND_SEED++;

                    struct timeval begin_utime;  // sec + usec.
                    struct timeval end_utime;  // sec + usec.

                    // insert each item into the bloom filter.
                    uint64_t num_inserts = 0;
                    gettimeofday(&begin_utime, NULL);
                    for(uint64_t k = 0; k < max_items; k++)
                    {
                        if (true == bloom_filter->Insert(item_array[k]))
                        {
                            num_inserts++;
                        }
                    }
                    gettimeofday(&end_utime, NULL);

                    bloom_filter->WriteFilterLog();

                    insert_result_txt << "false_positive_rate: " << FALSE_POSITIVE_RATE
                                      << ", table_occupancy_fraction: " << (1.0 * table_occupancy_fraction / 100)
                                      << ", max_inserts: " << max_items
                                      << ", num_ok_inserts: "  << num_inserts
                                      << ", insert_accesses: " << bloom_filter->CalculateInsertAccesses()
                                      << ", accesses_per_insert: " << (1.0 * bloom_filter->CalculateInsertAccesses() / num_inserts)
                                      << ", insert_time (s): " << ((end_utime.tv_sec - begin_utime.tv_sec) + ((end_utime.tv_usec - begin_utime.tv_usec) / 1000000.0))
                                      << ", insert_speed (mops): " << (1.0 * max_items / (1000000.0 * (end_utime.tv_sec - begin_utime.tv_sec) + (end_utime.tv_usec - begin_utime.tv_usec))) << endl;

                    cout << "false_positive_rate: " << FALSE_POSITIVE_RATE
                         << ", table_occupancy_fraction: " << (1.0 * table_occupancy_fraction / 100)
                         << ", max_inserts: " << max_items
                         << ", num_ok_inserts: "  << num_inserts
                         << ", insert_accesses: " << bloom_filter->CalculateInsertAccesses()
                         << ", accesses_per_insert: " << (1.0 * bloom_filter->CalculateInsertAccesses() / num_inserts)
                         << ", insert_time (s): " << ((end_utime.tv_sec - begin_utime.tv_sec) + ((end_utime.tv_usec - begin_utime.tv_usec) / 1000000.0))
                         << ", insert_speed (mops): " << (1.0 * max_items / (1000000.0 * (end_utime.tv_sec - begin_utime.tv_sec) + (end_utime.tv_usec - begin_utime.tv_usec))) << endl;

                    // lookup each item in the bloom filter.
                    uint32_t positive_lookup_percent = 0;
                    while(positive_lookup_percent <= 100)  // 0%, 50%, 100%.
                    {
                        uint64_t true_lookups = (uint64_t) (1.0 * max_tests * positive_lookup_percent / 100);
                        uint64_t false_lookups = max_tests - true_lookups;

                        // generate a test array.
                        uint64_t * test_array = new uint64_t [max_tests];
                        memset(test_array, 0, sizeof(uint64_t) * max_tests);
                        GenerateTestArray(item_array, max_items, test_array, max_tests, positive_lookup_percent);

                        bloom_filter->ResetMemoryAccesses();

                        // lookup each item in the bloom filter.
                        uint64_t num_lookups = 0;
                        gettimeofday(&begin_utime, NULL);
                        for(uint64_t k = 0; k < max_tests; k++)
                        {
                            if(true == bloom_filter->Lookup(test_array[k]))
                            {
                                num_lookups++;
                            }
                        }
                        gettimeofday(&end_utime, NULL);

                        lookup_result_txt << "false_positive_rate: " << FALSE_POSITIVE_RATE
                                          << ", max_lookups: " << max_tests
                                          << ", positive_lookup_percent: " << (1.0 * positive_lookup_percent / 100)
                                          << ", num_ok_lookups: " << num_lookups
                                          << ", num_nok_lookups: " << (max_tests - num_lookups)
                                          << ", lookup_accesses: " << bloom_filter->CalculateLookupAccesses()
                                          << ", accesses_per_lookup: " << (1.0 * bloom_filter->CalculateLookupAccesses() / max_tests)
                                          << ", false_positive_lookup_rate: " << ((false_lookups == 0) ? 0 : (1.0 * (num_lookups - true_lookups) / false_lookups))
                                          << ", lookup_time (s): " << ((end_utime.tv_sec - begin_utime.tv_sec) + ((end_utime.tv_usec - begin_utime.tv_usec) / 1000000.0))
                                          << ", lookup_speed (mops): " << (1.0 * max_tests / (1000000.0 * (end_utime.tv_sec - begin_utime.tv_sec) + (end_utime.tv_usec - begin_utime.tv_usec))) << endl;

                        cout << "false_positive_rate: " << FALSE_POSITIVE_RATE
                             << ", max_lookups: " << max_tests
                             << ", positive_lookup_percent: " << (1.0 * positive_lookup_percent / 100)
                             << ", num_ok_lookups: " << num_lookups
                             << ", num_nok_lookups: " << (max_tests - num_lookups)
                             << ", lookup_accesses: " << bloom_filter->CalculateLookupAccesses()
                             << ", accesses_per_lookup: " << (1.0 * bloom_filter->CalculateLookupAccesses() / max_tests)
                             << ", false_positive_lookup_rate: " << ((false_lookups == 0) ? 0 : (1.0 * (num_lookups - true_lookups) / false_lookups))
                             << ", lookup_time (s): " << ((end_utime.tv_sec - begin_utime.tv_sec) + ((end_utime.tv_usec - begin_utime.tv_usec) / 1000000.0))
                             << ", lookup_speed (mops): " << (1.0 * max_tests / (1000000.0 * (end_utime.tv_sec - begin_utime.tv_sec) + (end_utime.tv_usec - begin_utime.tv_usec))) << endl;

                        delete [] test_array;
                        positive_lookup_percent += 50;
                    }

                    bloom_filter->ResetMemoryAccesses();

                    // delete each item from the bloom filter.
                    uint32_t num_deletes = 0;
                    gettimeofday(&begin_utime, NULL);
                    for(size_t k = 0; k < max_items; k++)
                    {
                        if(true == bloom_filter->Delete(item_array[k]))
                        {
                            num_deletes++;
                        }
                    }
                    gettimeofday(&end_utime, NULL);

                    delete_result_txt << "false_positive_rate: " << FALSE_POSITIVE_RATE
                                      << ", table_occupancy_fraction: " << (1.0 * table_occupancy_fraction / 100)
                                      << ", max_deletes: " << max_items
                                      << ", num_ok_deletes: " << num_deletes
                                      << ", delete_accesses: " << bloom_filter->CalculateDeleteAccesses()
                                      << ", accesses_per_delete: " << (1.0 * bloom_filter->CalculateDeleteAccesses() / num_deletes)
                                      << ", delete_time (s): " << ((end_utime.tv_sec - begin_utime.tv_sec) + ((end_utime.tv_usec - begin_utime.tv_usec) / 1000000.0))
                                      << ", delete_speed (mops): " << (1.0 * max_items / (1000000.0 * (end_utime.tv_sec - begin_utime.tv_sec) + (end_utime.tv_usec - begin_utime.tv_usec))) << endl;

                    cout << "false_positive_rate: " << FALSE_POSITIVE_RATE
                         << ", table_occupancy_fraction: " << (1.0 * table_occupancy_fraction / 100)
                         << ", max_deletes: " << max_items
                         << ", num_ok_deletes: " << num_deletes
                         << ", delete_accesses: " << bloom_filter->CalculateDeleteAccesses()
                         << ", accesses_per_delete: " << (1.0 * bloom_filter->CalculateDeleteAccesses() / num_deletes)
                         << ", delete_time (s): " << ((end_utime.tv_sec - begin_utime.tv_sec) + ((end_utime.tv_usec - begin_utime.tv_usec) / 1000000.0))
                         << ", delete_speed (mops): " << (1.0 * max_items / (1000000.0 * (end_utime.tv_sec - begin_utime.tv_sec) + (end_utime.tv_usec - begin_utime.tv_usec))) << endl;

                    delete [] item_array;
                    delete bloom_filter;
                }

                table_occupancy_fraction += 10;
            }

            FALSE_POSITIVE_RATE /= 100;
        }

        SLOT_INDEX_BITS += 2;
    }

    insert_result_txt.close();
    lookup_result_txt.close();
    delete_result_txt.close();
}

int main()
{
    // test the bloom filter.
    TestBloomFilter();

    return 0;
}
