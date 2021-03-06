VERBOSITY_LEVEL = -1;
DataSet ds = ReadDataFile(temp.fasta);
DataSetFilter dsf = CreateFilter(ds,1,"","","");
HarvestFrequencies (vectorOfFrequencies, dsf,1,1,0);

global R = 1;
matrix = {{*,R*a,a,R*a}{R*a,*,R*a,a}{a,R*a,*,R*a}{R*a,a,R*a,*}};

Model HKY85Model = (matrix, vectorOfFrequencies, 1);

treeString = "(A,B,C);";
// Create tree
Tree fullTree = treeString;

LikelihoodFunction LF = (dsf, fullTree);
Optimize (res_LF,LF);
alt_LL = res_LF[1][0];

//fprintf(stdout,LF);

//constrain the length of one branch to be 0
fullTree.A.a :=0;

//reoptimize
Optimize (res_LF,LF);

//get new likelihood value
nullA_LL = res_LF[1][0];

//fprintf(stdout,LF);

//determine p value
lrtSTAT = 2*(alt_LL-nullA_LL);
pVal = 1-CChi2 (lrtSTAT, 1);

//print p-value to standard out
fprintf(stdout,pVal);
