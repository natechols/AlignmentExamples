
  $ printf ">seq1\nHEAGAWGHEE" > tmp1.fasta
  $ printf ">seq2\nPAWHEAE" > tmp2.fasta
  $ nwalign tmp1.fasta tmp2.fasta
  Sequence identity = 45.45%
  
  seq1  HEAGAWGHE-E
            ** ** *
  seq2  --P-AW-HEAE


  $ nwalign $TESTDIR/../data/hemoglobinA.fasta $TESTDIR/../data/hemoglobinB.fasta 
  Sequence identity = 43.24%
  
  hemoglobin_A  V-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DL
                * * *  *  * * ****     * * *** *     * *   *  * **
  hemoglobin_B  VHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDL
  
  
  hemoglobin_A  S-----HGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRV
                *      *   ** *****  *     ** *        ** **  ** *
  hemoglobin_B  STPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHV
  
  
  hemoglobin_A  DPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR
                ** ** **   *   ** *   **** * *   *  * *   *  ** 
  hemoglobin_B  DPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
