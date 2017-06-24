
#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <check.h>
#include "../src/nwalign.h"

START_TEST(test_nwalign_simple)
{
  char seq1[] = "HEAGAWGHEE";
  char seq2[] = "PAWHEAE";
  sequence_alignment *result = NeedlemanWunsch(seq1, seq2, "BLOSUM50", 8);
  ck_assert_msg(strcmp(result->seq1aligned, "HEAGAWGHE-E") == 0,
                "Expected HEAGAWGHE-E, got %s", result->seq1aligned);
  ck_assert_msg(strcmp(result->seq2aligned, "--P-AW-HEAE") == 0,
                "Expected --P-AW-HEAE, got %s", result->seq2aligned);
}
END_TEST

START_TEST(test_nwalign_hemoglobin)
{
#define HEM_A "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
#define HEM_B "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
#define HEM_A_ALN "V-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
#define HEM_B_ALN "VHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
  char seq1[] = HEM_A;
  char seq2[] = HEM_B;
  sequence_alignment *result = NeedlemanWunsch(seq1, seq2, "BLOSUM50", 8);
  ck_assert_msg(strcmp(result->seq1aligned, HEM_A_ALN) == 0,
                "Expected %s, got %s", HEM_A_ALN, result->seq1aligned);
  ck_assert_msg(strcmp(result->seq2aligned, HEM_B_ALN) == 0,
                "Expected %s, got %s", HEM_B_ALN, result->seq2aligned);
}
END_TEST

Suite *nwalign_suite(void) {
  Suite *s;
  TCase *tc_1, *tc_2;

  s = suite_create("nwalign");
  tc_1 = tcase_create("simple");
  tcase_add_test(tc_1, test_nwalign_simple);
  tc_2 = tcase_create("hemoglobin");
  tcase_add_test(tc_2, test_nwalign_hemoglobin);
  suite_add_tcase(s, tc_1);
  suite_add_tcase(s, tc_2);
  return s;
}

int main(void)
{
  int number_failed;
  Suite *s;
  SRunner *sr;
  s = nwalign_suite();
  sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
