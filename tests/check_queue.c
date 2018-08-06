#include "check_all.h"
#include "../src/queue.h"


static KC__MemAllocator* ma;
static KC__Queue* q;

static void setup() {
    ma = KC__mem_allocator_create(1000000);
    q = KC__queue_create(ma, 3);
}

static void teardown() {
    KC__queue_free(ma, q);
    KC__mem_allocator_free(ma);
}

START_TEST(test_function)
    {
        ck_assert(KC__queue_capacity(q) == 3);
        ck_assert(KC__queue_length(q) == 0);
        ck_assert(KC__queue_is_empty(q));

        int d[3];

        // Enqueue and dequeue _i times to ensure circular works.
        for (int i = 0; i < _i; i++) {
            bool r = KC__queue_enqueue(q, &(d[0]));
            ck_assert(r);
            int* v = KC__queue_dequeue(q);
            ck_assert(v == &(d[0]));
            ck_assert(KC__queue_is_empty(q));
        }

        // Enqueue 3 items, the front should be the address of d[0].
        for (int i = 0; i < 3; i++) {
            bool r = KC__queue_enqueue(q, &(d[i]));
            ck_assert(r);
            ck_assert(KC__queue_length(q) == i + 1);
            ck_assert(KC__queue_front(q) == &(d[0]));
        }

        ck_assert(KC__queue_length(q) == 3);
        ck_assert(KC__queue_is_full(q));

        bool r = KC__queue_enqueue(q, &(d[0]));
        ck_assert_msg(!r, "Queue is full, enqueue should fail.");

        // Dequeue 3 items.
        for (int i = 0; i < 3; i++) {
            ck_assert(KC__queue_front(q) == &(d[i]));
            int* v = KC__queue_dequeue(q);
            ck_assert(v == &(d[i]));
        }

        ck_assert(KC__queue_length(q) == 0);
        ck_assert(KC__queue_is_empty(q));

        int* v = KC__queue_dequeue(q);
        ck_assert_msg(v == NULL, "Queue is empty, dequeue should return NULL");
    }
END_TEST

Suite* queue_suite() {
    TCase* tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core, setup, teardown);
    tcase_add_loop_test(tc_core, test_function, 0, 10);

    Suite* s = suite_create("Queue");
    suite_add_tcase(s, tc_core);

    return s;
}