QUnit.test("Simple alignment test", function (assert) {
  var m = scoreMatrices["BLOSUM50"];
  var aln = nwalign.NeedlemanWunsch("HEAGAWGHEE", "PAWHEAE", m, 8);
  console.log(aln.seq1aligned);
  console.log(aln.seq2aligned);
  assert.ok(aln.seq1aligned == "HEAGAWGHE-E", "SEQ1 Passed");
  assert.ok(aln.seq2aligned == "--P-AW-HEAE", "SEQ2 passed");
  assert.ok(aln.identity() - 0.4545 < 0.0001, "Identity passed");
  console.log(aln.format("seq1", "seq2"));
});

QUnit.test("Hemoglobin alignment", function (assert) {
  var m = scoreMatrices["BLOSUM50"];
  var aln = nwalign.NeedlemanWunsch("VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR", "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH", m, 8);
  assert.ok(aln.seq1aligned == "V-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR", "Passed");
  assert.ok(aln.seq2aligned == "VHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH", "Passed");
  console.log(aln.identity());
  assert.ok(aln.identity() - 0.4324 < 0.0001, "Identity passed");
});
