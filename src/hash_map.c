/*
 * This file is part of CHTKC.
 *
 * CHTKC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CHTKC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CHTKC.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: Jianan Wang
 */


#include <pthread.h>

#ifdef __APPLE__
#include "pthread_barrier.h"
#endif // __APPLE__

#include <stdlib.h>
#include <string.h>
#include "hash_map.h"
#include "logging.h"
#include "assert.h"
#include "utils.h"


typedef struct {
    KC__node_id_t next;

    /** The count equaling to 0 means this node has not added into hash table. */
    KC__count_t count;

    KC__unit_t kmer[];
} KC__HashMapNode;


typedef struct {
    KC__node_id_t start_id;
    KC__node_id_t end_id;

    /**
     * The current node is pre-fetched, so that there always be an available node when adding a K-mer unless the keys
     * of hash map have been locked. If the K-mer already exists, current_id is still valid, else current_id should be
     * updated, if no new node can be found, the hash map should lock keys.
     */
    KC__node_id_t current_id;

    /**
     * The field is used to get a new node. Only nodes between start_id and next_id have already added into hash table.
     * If a thread has already used out all nodes in its block, it can request nodes from other blocks by accessing
     * next_id field.
     */
    KC__node_id_t next_id;

    bool synced;
} KC__HashMapNodeBlock;


struct KC__HashMap {
    KC__node_id_t* table;
    size_t table_capacity;

    /** The node at position 0 is reserved as NULL. */
    KC__HashMapNode* nodes;
    size_t node_size;
    size_t kmer_size;
    size_t kmer_width;

    KC__HashMapNodeBlock** blocks;
    size_t blocks_count;

    bool keys_locked;
    pthread_barrier_t barrier;
};


KC__HashMap* KC__hash_map_create(KC__MemAllocator* ma, size_t K, size_t threads_count) {
    KC__HashMap* hm = (KC__HashMap*)KC__mem_alloc(ma, sizeof(struct KC__HashMap), "hash map");

    hm->blocks_count = threads_count;
    hm->blocks = (KC__HashMapNodeBlock**)KC__mem_alloc(ma, sizeof(KC__HashMapNodeBlock*) * hm->blocks_count, "hash map blocks array");
    for (size_t i = 0; i < hm->blocks_count; i++) {
        hm->blocks[i] = (KC__HashMapNodeBlock*)KC__mem_aligned_alloc(ma, sizeof(KC__HashMapNodeBlock), "hash map block");
    }

    hm->kmer_width = KC__calculate_kmer_width(K);
    hm->kmer_size = KC__calculate_kmer_size(K);
    hm->node_size = sizeof(KC__HashMapNode) + hm->kmer_size;

    const size_t mem_limit = KC__mem_available(ma);

    size_t nodes_count_limit = mem_limit / (hm->node_size * 3 + sizeof(KC__node_id_t) * 4) * 3;
    if (nodes_count_limit > KC__NODE_ID_MAX) {
        LOGGING_WARNING("The count of nodes to be allocated is too large: %zu.", nodes_count_limit);
        nodes_count_limit = KC__NODE_ID_MAX;
        LOGGING_WARNING("Reduce the count of nodes to %zu.", nodes_count_limit);
    }
    KC__node_id_t nodes_count = (KC__node_id_t)nodes_count_limit;
    const size_t nodes_mem = hm->node_size * nodes_count;

    const size_t table_mem_limit = mem_limit - nodes_mem;
    const size_t table_capacity_limit = table_mem_limit / sizeof(KC__node_id_t);
    hm->table_capacity = KC__max_prime_number(table_capacity_limit);
    const size_t table_mem = sizeof(KC__node_id_t) * hm->table_capacity;
    hm->table = (KC__node_id_t*)KC__mem_aligned_alloc(ma, table_mem, "hash map table");

    hm->nodes = (KC__HashMapNode*)KC__mem_aligned_alloc(ma, nodes_mem, "hash map nodes");
    const KC__node_id_t step = nodes_count / (KC__node_id_t)(hm->blocks_count);
    for (size_t i = 0; i < hm->blocks_count; i++) {
        KC__HashMapNodeBlock* block = hm->blocks[i];
        block->start_id = 1 + step * (KC__node_id_t)i;
        block->end_id = (i == hm->blocks_count - 1) ? (nodes_count) : (1 + step * (KC__node_id_t)(i + 1));
    }

    LOGGING_DEBUG("        Hash table capacity: %zu (limit: %zu)", hm->table_capacity, table_capacity_limit);
    LOGGING_DEBUG("          Hash table memory: %zu", table_mem);
    LOGGING_DEBUG("                Nodes count: %zu", nodes_count);
    LOGGING_DEBUG("               Nodes memory: %zu", nodes_mem);
    LOGGING_DEBUG("Hash table and nodes memory: %zu (limit: %zu)", table_mem + nodes_mem, mem_limit);
    for (size_t i = 0; i < hm->blocks_count; i++) {
        KC__HashMapNodeBlock* block = hm->blocks[i];
        LOGGING_DEBUG("Nodes block #%zu (start: %zu, end: %zu, length: %zu)", i, block->start_id, block->end_id, block->end_id - block->start_id);
    }

    pthread_barrier_init(&(hm->barrier), NULL, (unsigned int)threads_count);

    KC__hash_map_clear(hm);

    return hm;
}

void KC__hash_map_free(KC__MemAllocator* ma, KC__HashMap* hm) {
    for (size_t i = 0; i < hm->blocks_count; i++) {
        KC__mem_free(ma, hm->blocks[i]);
    }
    KC__mem_free(ma, hm->blocks);

    KC__mem_free(ma, hm->nodes);
    KC__mem_free(ma, hm->table);

    pthread_barrier_destroy(&(hm->barrier));

    KC__mem_free(ma, hm);
}

size_t KC__hash_map_max_key_count(const KC__HashMap* hm) {
    // The node at position 0 is reserved as NULL, which cannot hold a valid key.
    return hm->blocks[hm->blocks_count - 1]->end_id - 1;
}

void KC__hash_map_set_table_capacity(KC__HashMap* hm, size_t capacity) {
    LOGGING_WARNING("Set table capacity to %zu (should only be used for tests)", capacity);
    KC__ASSERT(capacity <= hm->table_capacity);
    hm->table_capacity = capacity;
}

void KC__hash_map_lock_keys(KC__HashMap* hm) {
    LOGGING_WARNING("Set hash table key locked (should only be used for tests)");
    hm->keys_locked = true;
}

static inline bool KC__hash_map_kmers_equal(const KC__HashMap* hm, const KC__unit_t* kmer_1, const KC__unit_t* kmer_2) {
    for (size_t i = 0; i < hm->kmer_width; i++) {
        if (kmer_1[i] != kmer_2[i]) {
            return false;
        }
    }
    return true;
}

static inline void KC__hash_map_copy_kmer(KC__HashMap* hm, KC__unit_t* dest, const KC__unit_t* src) {
    memcpy(dest, src, hm->kmer_size);
}

static inline KC__HashMapNode* KC__hash_map_get_node(const KC__HashMap* hm, KC__node_id_t node_id) {
    char* node = (char*)(hm->nodes);
    node += hm->node_size * (size_t)node_id;
    return (KC__HashMapNode*)node;
}

static inline KC__node_id_t KC__hash_map_request_node(KC__HashMap* hm, size_t n) {
    KC__HashMapNodeBlock* block = hm->blocks[n];
    KC__node_id_t node_id;

    do {
        node_id = block->next_id;
        if (node_id == block->end_id) {
            return KC__NODE_ID_NULL;
        }
    } while (!__sync_bool_compare_and_swap(&(block->next_id), node_id, node_id + 1));

    KC__HashMapNode* node = KC__hash_map_get_node(hm, node_id);
    node->count = 0;

    return node_id;
}

static inline KC__node_id_t KC__hash_map_polling_request_node(KC__HashMap* hm, size_t n) {
    KC__node_id_t node_id = KC__hash_map_request_node(hm, n);
    if (node_id != KC__NODE_ID_NULL) {
        return node_id;
    }

    size_t m = hm->blocks_count;
    for (size_t i = 0; i < m - 1; i++) {
        n++;
        if (n == m)
            n = 0;

        node_id = KC__hash_map_request_node(hm, n);
        if (node_id != KC__NODE_ID_NULL) {
            return node_id;
        }
    }
    return KC__NODE_ID_NULL;
}

typedef struct {
    KC__HashMap* hash_map;
    size_t n;
    size_t start;
    size_t end;
} KC__HashMapClearTableParam;

static void* KC__hash_map_clear_table(void* ptr) {
    KC__HashMapClearTableParam* param = (KC__HashMapClearTableParam*)ptr;
    KC__HashMap* hm = param->hash_map;
    size_t n = param->n;
    size_t start = param->start;
    size_t end = param->end;
    LOGGING_DEBUG("Hash table clear #%zu from %zu to %zu (length: %zu)", n, start, end, end - start);
    for (size_t i = param->start; i < param->end; i++) {
        hm->table[i] = KC__NODE_ID_NULL;
    }
    pthread_exit(NULL);
}

void KC__hash_map_clear(KC__HashMap* hm) {
    hm->keys_locked = false;

    for (size_t i = 0; i < hm->blocks_count; i++) {
        KC__HashMapNodeBlock* block = hm->blocks[i];
        block->next_id = block->start_id;
        block->current_id = KC__NODE_ID_NULL;
        block->synced = false;
    }

    // The count of threads used to clear hash table equals to blocks count.
    size_t clear_table_threads_count = hm->blocks_count;
    pthread_t threads[clear_table_threads_count];
    KC__HashMapClearTableParam params[clear_table_threads_count];
    size_t step = hm->table_capacity / clear_table_threads_count;
    for (size_t i = 0; i < clear_table_threads_count; i++) {
        size_t start = i * step;
        size_t end = (i == clear_table_threads_count - 1) ? (hm->table_capacity) : ((i + 1) * step);
        params[i].hash_map = hm;
        params[i].n = i;
        params[i].start = start;
        params[i].end = end;
        pthread_create(&(threads[i]), NULL, KC__hash_map_clear_table, &(params[i]));
    }

    for (size_t i = 0; i < clear_table_threads_count; i++) {
        pthread_join(threads[i], NULL);
    }
}

static inline size_t KC__hash_map_hash_function(const KC__HashMap* hm, const KC__unit_t* kmer) {
    size_t n = 0;
    for (size_t i = 0; i < hm->kmer_width; i++) {
        n += kmer[i];
    }
    return n % hm->table_capacity;
}

/**
 * Add K-mer to collision list (may be part of the list) specified by pointer to a node id.
 * @param hm The hash map.
 * @param kmer The K-mer to be added.
 * @param list Specify the head of the (sub-) collision list, will be updated before return.
 * @return If the K-mer already exists in the collision list, the node id will be returned, and list will be updated to
 * a pointer to this node, else the tail (KC__NODE_ID_NULL) of collision list will be returned and list will be updated
 * to the corresponding pointer.
 */
static inline KC__node_id_t KC__hash_map_collision_list_add_kmer(KC__HashMap* hm, KC__node_id_t** list, const KC__unit_t* kmer) {
    KC__node_id_t node_id;
    KC__node_id_t* p = *list;
    while (true) {
        node_id = *p;

        if (node_id == KC__NODE_ID_NULL) {
            break;
        }

        KC__HashMapNode* node = KC__hash_map_get_node(hm, node_id);
        if (KC__hash_map_kmers_equal(hm, node->kmer, kmer)) {
            KC__count_t count;
            do {
                count = node->count;
                if (count == KC__COUNT_MAX)
                    break;
            } while (!__sync_bool_compare_and_swap(&(node->count), count, count + 1));

            break;
        }
        p = &(node->next);
    }

    *list = p;
    return node_id;
}

bool KC__hash_map_add_kmer(KC__HashMap* hm, size_t n, const KC__unit_t* kmer) {
    KC__HashMapNodeBlock* block = hm->blocks[n];

    if (!(block->synced) && (block->current_id == KC__NODE_ID_NULL)) {
        block->current_id = KC__hash_map_polling_request_node(hm, n);
        if (block->current_id == KC__NODE_ID_NULL) {
            hm->keys_locked = true;
            LOGGING_DEBUG("Set hash map keys locked.");
        }
    }

    // For multi-threaded situation, when the keys are locked by one of the working threads, the other threads may
    // not yet knew the change, so they need to sync once. As the failure of adding K-mer may occur a long time later
    // (if one thread adds many existed K-mers continuously), checking if the keys of the hash map are locked is a
    // good signal.

    if (!(block->synced) && hm->keys_locked) {
        pthread_barrier_wait(&(hm->barrier));
        block->synced = true;
        LOGGING_DEBUG("Block #%zu synced (keys locked).", n);
    }

    size_t table_idx = KC__hash_map_hash_function(hm, kmer);

    KC__node_id_t* collision_list = &(hm->table[table_idx]);
    KC__node_id_t node_id = KC__hash_map_collision_list_add_kmer(hm, &collision_list, kmer);

    if (node_id != KC__NODE_ID_NULL) {
        return true;
    }

    // If some thread has set keys_locked, assume this thread noticed the change here and return false, while the other
    // thread has not seen the change and is adding a new node to hash table, then it may cause inconsistency.
    // Checking if this thread has been synced is important.
    if (block->synced && hm->keys_locked) {
        return false;
    }

    KC__HashMapNode *node = KC__hash_map_get_node(hm, block->current_id);
    KC__hash_map_copy_kmer(hm, node->kmer, kmer);
    node->count = 1;
    node->next = KC__NODE_ID_NULL;

    do {
        node_id = KC__hash_map_collision_list_add_kmer(hm, &collision_list, kmer);
        if (node_id != KC__NODE_ID_NULL) {
            // Mark the node invalid.
            node->count = 0;
            return true;
        }
    } while (!__sync_bool_compare_and_swap(collision_list, node_id, block->current_id));

    block->current_id = KC__NODE_ID_NULL;

    return true;
}

void KC__hash_map_finish_adding_kmers(KC__HashMap* hm, size_t n) {
    KC__HashMapNodeBlock* block = hm->blocks[n];

    if (!(block->synced)) {
        pthread_barrier_wait(&(hm->barrier));
        block->synced = true;
        LOGGING_DEBUG("Block #%zu synced (adding finished).", n);
    }
}

void KC__hash_map_export(KC__HashMap* hm, size_t n, KC__HashMapExportCallback callback, void* data, size_t* exported_count) {
    size_t ec = 0;

    KC__ASSERT(n < hm->blocks_count);

    KC__HashMapNodeBlock* block = hm->blocks[n];
    for (KC__node_id_t i = block->start_id; i < block->next_id; i++) {
        KC__HashMapNode* node = KC__hash_map_get_node(hm, i);
        if (node->count != 0) {
            callback(node->kmer, node->count, data);
            ec++;
        } else {
            LOGGING_DEBUG("Block #%zu (%zu-%zu) node id: %zu count equals to 0.", i, block->start_id, block->end_id, i);
        }
    }

    if (exported_count != NULL) {
        *exported_count = ec;
    }
}