 /* This program uses King et al.'s matrix inversion technique to calculate the aerosol number density.
	It requires the Mie scattering coefficients calculated by bhmie.c for each wavelength used
	over the radius of interest.
    This program assumes there is spectral data from 6 wavelengths: 368nm, 500nm, 673nm, 778nm, 870nm, 940nm.
    The radius of aerosols affecting each wavelength are determined from Mie theory: Qsca >= 3
    This gives (wavelength in nm, radii in micrometers):
    
    Wavelength	r_min	r_max	r_ave	length of range	index,i
    ===========================================================
    368		0.174	0.310	0.247		        147		            0
    500		0.236	0.412	0.324			177			    1
    673		0.318	0.567	0.443			250			    2
    778		0.367	0.655	0.511			289			    3
    870		0.411	0.773	0.572			323			    4
    940		0.443	0.792	0.618			350	    	            5
    
    These overlapping radius ranges determine the structure of King's A matrix (eq. 5).
    Excluding the water vapor absorption channels, form subintervals:
    
    index,j	radius range	include in g_i for i = 		r_ave
    =================================================================
    0		0.174-0.236	   0			        0.205
    1		0.236-0.310	   0  1	                        0.273
    2		0.310-0.318	      1                         0.314
    3		0.318-0.367	      1  2  		        0.343
    4		0.367-0.411	      1  2  3		        0.389
    5		0.411-0.412	      1  2  3  4	        0.412
    6		0.412-0.443	         2  3  4	        0.428
    7		0.443-0.567	         2  3  4  5	        0.505
    8		0.567-0.655		    3  4  5	        0.611
    9		0.655-0.733	               4  5  	        0.694
    10		0.733-0.792		          5             0.763
   
    Now the index i represents the wavelength (row) and the index j represents the radius subinterval (column)
    
	         | A_00  A_01    0     0     0     0     0     0     0       0       0         |
	         |   0   A_11  A_12  A_13  A_14  A_15    0     0     0       0       0         |
      A_ij = 	 |   0     0     0   A_23  A_24  A_25  A_26   A_27   0       0       0         |
	         |   0     0     0     0   A_34  A_35  A_36   A_37   A_38    0       0         |
	         |   0     0     0     0     0   A_45  A_46   A_47   A_48   A_49     0         |
           	 |   0	   0     0     0     0     0     0    A_57   A_58   A_59    A_510      |
*/

// c++ includes
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <new>
#include <cmath>
#include <vector>

using namespace std;	// required for c++ file routines

#define N0 1.0E12	// an initial number estimate guess
#define PI 3.1415192654		// the constant PI
#define CONV_FACT 1.0E-4	/* convert microns to cm */

// required definition of matrix at top of source code
typedef vector<vector<double> > matrix;

// function prototype definitions to be placed before main() source code
void Sinvmtrx(matrix &a, matrix &b);	// inverts a matrix - assumes non-zero diagonal elements: b = a^-1
void mtrxmult(matrix &X,matrix &Y,matrix &Z);	// multiplies matrices: Z = X * Y
void Tmtrxmult(matrix &X,matrix &Y,matrix &Z);	// multiplies the transpose of a matrix with a matrix: Z = X^T * Y
void mtrxprnt(matrix &a);	// prints a matrix to the screen - for troubleshooting
// Function Adepth calculates the optical depth over sub-range of detector as defined in the A matrix
double Adepth(double rad[], double Qext[], double npart[], int length, int offset, double step);
// Function optdepth calculates the optical depth over full range of detector
double optdepth(double rad[], double Qext[], double npart[], int length, int offset, double step);
// Function initQ initializes the mie scattering extinction vectors
void initQ(double Q[], int length, string fname);
// Function initn initializes the number density vector
void initn(double r[], double n[]);
// Function initH intializes the smoothing matrix H
void initH(matrix &H);
// Function nadjust modifies the number density vector n(r)=n(r)f(r)
void nadjust(double n[], int len, int offset, double f);
// Function optimum_gamma finds the minimum value of gamma for which all elements of solution vector f are positive
double optimum_gamma(matrix &AtCiA, matrix &H, matrix &AtCi, double aod[]);

int main()
{
	string mystr,fname;
	char *fname1,*fname2,*fname4,*fname5,dataline[500];
	int i,j,k,l,channel;
	double day;
	int wave[6]={368,500,673,778,870,840};	// detector center wavelengths
	double gamma=1.0E-3;		// Lagrangian multiplier
	double r[619],n[619];		// Pointers for radius and number density arrays
	double aodm[6], eaodm[6], aodc[6];	// measured and calculated aod at 6 wavelengths
	double nest[6];		// number estimate at each of the 6 wavelengths
	fstream aodfile,outfile,narray,summary;
	double h=0.001*CONV_FACT;		// step size for numerical integration
	// holds Qext from mie calculations
	double Q368[137],Q500[177],Q673[250],Q778[289],Q870[323],Q940[350];
	matrix H(11,vector<double>(11,0.0));		// smoothing matrix minimizes second derivatives
	matrix A(6,vector<double>(11,0.0));		// optical depth matrix
	matrix C(6,vector<double>(6,0.0));		// diagonal weighting matrix of measurement uncertainties
	matrix Cinv(6,vector<double>(6,0.0));		// the inverse of C
	matrix B(11,vector<double>(11,0.0));		// B = A^T * C^-1 * A + gamma * H
	matrix Binv(11,vector<double>(11,0.0));		// the inverse of B
	matrix AtCinv(11,vector<double>(6,0.0));		// A transposed times inverse of C
	matrix AtCinvA(11,vector<double>(11,0.0));	// A transposed times inverse of C times A
	matrix BinvAtCinv(11,vector<double>(6,0.0));	// inverse of B times A transposed times inverse of C times A
	double f[11];					// solution vector f
	for(i=0;i<6;i++){aodc[i]=0.0;}			// initialize calculated aerosol optical depths to zero


	// fill detector arrays in this order: 368nm, 500nm, 673nm, 778nm, 870nm, 940nm
	// arrays for radius range parameters
	double ri[6]={0.242,0.324,0.443,0.511,0.572,0.618};	// mean radius for each detector's radius range (microns)
	int roff[6]={0,62,144,193,237,269};			// offsets from smallest radius to start of each detector's radius range
	int len[6]={137,177,250,289,323,350};	// length of radius ranges for 368nm, 500nm, 673nm, 778nm, 870nm, 940nm detectors
	// arrays for radius interval parameters
	double rbar[11]={0.205,0.273,0.314,0.343,0.389,0.412,0.428,0.505,0.611,0.694,0.763};	// mean radius for each radii sub-interval (microns)
	int length[11]={63,75,9,50,45,2,32,125,89,79,69};			// length of radii sub-intervals
	int roffset[11]={0,62,136,144,193,237,238,269,393,481,559};	// offset for radii sub-intervals

	// initialize mie extinction vectors
	fname="mie368.csv";
	initQ(Q368,len[0],fname);
	fname="mie500.csv";
	initQ(Q500,len[1],fname);
	fname="mie673.csv";
	initQ(Q673,len[2],fname);
	fname="mie778.csv";
	initQ(Q778,len[3],fname);
	fname="mie870.csv";
	initQ(Q870,len[4],fname);
	fname="mie940.csv";
	initQ(Q940,len[5],fname);
	

	cout << "\nThis program uses King's method to fit an aerosol size distribution to 6 channels of a spectrometer.\n";

	// Initialize the radius (microns) and number density matrices
	for(i=0;i<619;i++){r[i]=0.221 + i/1000.0;}
	initn(r,n);
	// Initialize the smoothing matrix H
	initH(H);
	cout << "Enter measured AOD file to process (e.g., AOD.csv): ";
	getline(cin,mystr);
	fname1 = new char [mystr.size()+1];	/* c_string method */
	strcpy(fname1,mystr.c_str());

	cout << "Enter output file name for aerosol number densities (e.g., 2010aer-number.csv): ";
	getline(cin,mystr);
	fname2 = new char [mystr.size()+1];	/* c_string method */
	strcpy(fname2,mystr.c_str());

	// cout << "Enter output file name for 930/940nm aerosol optical depths (e.g., 2010aod-930-940.csv): ";
	//getline(cin,mystr);
	//fname3 = new char [mystr.size()+1];	/* c_string method */
	// strcpy(fname3,mystr.c_str());
     
	cout << "Enter output file name for aerosol number array (e.g., 2010num-array.csv): ";
	getline(cin,mystr);
	fname4 = new char [mystr.size()+1];	/* c_string method */
	strcpy(fname4,mystr.c_str());

	cout << "Enter output file name for summary information(e.g., 2010summary.txt): ";
	getline(cin,mystr);
	fname5 = new char [mystr.size()+1];	/* c_string method */
	strcpy(fname5,mystr.c_str());

	aodfile.open(fname1);	// could use myfile.open(fname); 
	if (!aodfile.is_open()){
		cout << "\nUnable to open file " << fname1 << "\n";
		return 0;
	}
	outfile.open(fname2,ios::out);	// could use myfile.open(fname); 
	if (!outfile.is_open()){
		cout << "\nUnable to open file " << fname2 << "\n";
		return 0;
	}
	
	//h2ofile.open(fname3,ios::out);	// could use myfile.open(fname); 
	//if (!h2ofile.is_open()){
	//	cout << "\nUnable to open file " << fname3 << "\n";
	//	return 0;
	//}
	
	narray.open(fname4,ios::out);	// could use myfile.open(fname); 
	if (!narray.is_open()){
		cout << "\nUnable to open file " << fname4 << "\n";
		return 0;
	}
	summary.open(fname5,ios::out);	// could use myfile.open(fname); 
	if (!summary.is_open()){
		cout << "\nUnable to open file " << fname5 << "\n";
		return 0;
	}
	outfile << "Day";		// write header line to output file - day
	for(i=0;i<6;i++){outfile << "," << ri[i];} // now write mean radiii 
	outfile << endl;
	//h2ofile << "Day,930aod,940aod" << endl;		// write header line to h2o aod file
	narray << "r,n" << endl; // write header line to number density file
	
	i=0;
	while(!aodfile.eof())	// use this line for final runs                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//while(i<1)	// this line for testing - do only one set of data for testing
	{
		// reset the number density matrix for each day
		initn(r,n);
		if(i==0){getline(aodfile,mystr);}	// toss out header line of aod data file
		getline(aodfile,mystr);			// get a line of data
		strcpy(dataline,mystr.c_str());		// extract data fields
		sscanf(dataline,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&day,&aodm[0],&aodm[1],&aodm[2],&aodm[3],&aodm[4],&aodm[5]);
		cout << day << endl;
		for(j=0;j<6;j++){
			eaodm[j]=0.02*aodm[j];
			C[j][j]=eaodm[j];} // 368nm, 500nm, 673nm, 778nm, 870nm, 940nm measurement uncertainties
		Sinvmtrx(C,Cinv);
		k=0;	// will do k iterations for each set of measurements
		do {
			// Calculate elements of the optical depth matrix A for each iteration
			// 368 nm
			for(j=0;j<2;j++){A[0][j]=Adepth(r,Q368,n,length[j],roffset[j],h);}
			// 500 nm
			for(j=1;j<6;j++){A[1][j]=Adepth(r,Q500,n,length[j],roffset[j],h);}
            		// 673 nm
			for(j=3;j<8;j++){A[2][j]=Adepth(r,Q673,n,length[j],roffset[j],h);}
			// 778 nm
			for(j=4;j<9;j++){A[3][j]=Adepth(r,Q778,n,length[j],roffset[j],h);}
			// 870 nm
			for(j=5;j<10;j++){A[4][j]=Adepth(r,Q870,n,length[j],roffset[j],h);}
			// 862 nm
			for(j=7;j<11;j++){A[5][j]=Adepth(r,Q940,n,length[j],roffset[j],h);}
			// Calculate A^T * C^-1
			Tmtrxmult(A,Cinv,AtCinv);
			// Calculate A^T * C^-1 * A
			mtrxmult(AtCinv,A,AtCinvA);
			gamma=optimum_gamma(AtCinvA,H,AtCinv,aodm);
			// Calculate B
			for(l=0;l<11;l++){for(j=0;j<11;j++){B[l][j]=AtCinvA[l][j]+gamma*H[l][j];}}
			// Calculate B^-1
			Sinvmtrx(B,Binv);
			// Calculate B^-1 * A^T * C^-1
			mtrxmult(Binv,AtCinv,BinvAtCinv);
			// Calculate corrections f[9]
			for(j=0;j<11;j++){f[j]=0.0;}  // first reset f vector
			for(j=0;j<11;j++){
				for(l=0;l<6;l++){
					f[j]+=BinvAtCinv[j][l]*aodm[l];
				}
				// and update n(r) intervals
				nadjust(n,length[j],roffset[j],f[j]);
			}
			k++;
// end loop here for final runs, otherwise end below for testing
		} while((k<20));

		// Now calculate optical depths for wavelength ranges                             
		// Do 368 nm, channel 0, len=137, range of i: 0-137, offset=0
		channel=0;
		aodc[channel]=optdepth(r,Q368,n,len[channel],roff[channel],h);
		
		// Do 500 nm, channel 1, len=177, range of i: 62-238, offset=62
		channel=1;
		aodc[channel]=optdepth(r,Q500,n,len[channel],roff[channel],h);
		
		// Do 673 nm, channel 2, len=250, range of i: 144-393, offset=144
		channel=2;
		aodc[channel]=optdepth(r,Q673,n,len[channel],roff[channel],h);
		
		// Do 778 nm, channel 3, len=289, range of i: 193-481, offset=193                  
		channel=3;
		aodc[channel]=optdepth(r,Q778,n,len[channel],roff[channel],h);
		
		// Do 870 nm, channel 4, len=323, range of i: 237-559, offset=237
		channel=4;
		aodc[channel]=optdepth(r,Q870,n,len[channel],roff[channel],h);
		
		// Do 940 nm, channel 5, len=307, range of i: 269-618, offset=269
		channel=5;
		aodc[channel]=optdepth(r,Q940,n,len[channel],roff[channel],h);
		

		nest[0]=n[68];	// number estimate for 368nm at average radius
		nest[1]=n[150];	// number estimate for 500nm at average radius
		nest[2]=n[269];	// number estimate for 673nm at average radius
		nest[3]=n[337];	// number estimate for 778nm at average radius
		nest[4]=n[398];	// number estimate for 870nm at average radius
		nest[5]=n[444];	// number estimate for 940nm at average radius

		cout << "\nSummary results for day " << day << " (" << k << " iterations):\n" << "Wavelength\tAOD (measured)\tUncertainty\tAOD (calculated)\tr (microns)\tn(r)\n";
		for(i=0;i<6;i++){
		    cout << wave[i] << "\t\t" << aodm[i] << "\t" << eaodm[i] << "\t" << aodc[i] << "\t\t" << ri[i] << "\t\t" << nest[i] << endl;
		}
		summary << "\nSummary results for day " << day << " (" << k << " iterations):\n" << "Wavelength\tAOD (measured)\tUncertainty\tAOD (calculated)\tr (microns)\tn(r)\n";
		for(i=0;i<6;i++){
		    summary << wave[i] << "\t\t" << aodm[i] << "\t" << eaodm[i] << "\t" << aodc[i] << "\t\t" << ri[i] << "\t\t" << nest[i] << endl;
		}
// end loop here for testing, otherwise end it on line above
//		} while((k<20));

		if(!aodfile.eof()){
		    outfile << day;                                                           
		    for(j=0;j<6;j++){outfile << "," << nest[j];}                        
		    outfile << endl;
			//h2ofile << day << "," << aodc[5] << "," << aodc[6] << endl;
		}
//		i++;	// for testing purposes only
	}
	for(i=0;i<619;i++){narray << r[i] << "," << n[i] << endl;}
	summary.close();
	narray.close();
	aodfile.close();
	outfile.close();
	//h2ofile.close();
	cout << "Processed aod file " << fname1 << "\n";
	cout << "Output particle data written to " << fname2 << "\n";
	//cout << "Output aod data for 930/940nm written to " << fname3 << endl;
	cout << "Final day's n(r) data written to " << fname4 << endl;
	cout << "Summary information written to " << fname5 << endl << endl;

	return 0;
}

// functions to be placed at end of main() source code
double Adepth(double rad[], double Qext[], double npart[], int length, int offset, double step)
/* Calculates the optical depth over sub-range of detector as defined in the A matrix */
{
	int i,j;
	double f=0.0;
	for (i=0; i<length; i++){
		j=i+offset;
		f=f+PI*pow(rad[j]*CONV_FACT,2)*Qext[j]*npart[j];
	}
	f=f-0.5*(PI*pow(rad[offset]*CONV_FACT,2)*Qext[offset]*npart[offset]+PI*pow(rad[length+offset-1]*CONV_FACT,2)*Qext[length+offset-1]*npart[length+offset-1]);
	f*=step;
	return f;
}
double optdepth(double rad[], double Qext[], double npart[], int length, int offset, double step)
/* Calculates the optical depth over full range of detector */
{
	int i,j;
	double f=0.0;
	for (i=0; i<length; i++){
		j=i+offset;
		f=f+PI*pow(rad[j]*CONV_FACT,2)*Qext[i]*npart[j];
	}
	f=f-0.5*(PI*pow(rad[offset]*CONV_FACT,2)*Qext[0]*npart[offset]+PI*pow(rad[length+offset-1]*CONV_FACT,2)*Qext[length-1]*npart[length+offset-1]);
	f*=step;
	return f;
}
void initQ(double Q[], int length, string fname)
// extracts the relevant Mie scattering values from the Mie files generated by the bhmie.c program
{
	int i;
	double r, Qsca, Qabs, Gsca;	// data in Mie files but not needed
	char *fname1,dataline[100];
	string mystr;
	fstream miefile;
	fname1 = new char [fname.size()+1];	/* c_string method */
	strcpy(fname1,fname.c_str());

	miefile.open(fname1);	// could use myfile.open(fname); 
	if (miefile.is_open())
	{
		i=0;
		while(i<length)
		{
			getline(miefile,mystr);
			strcpy(dataline,mystr.c_str());
			sscanf(dataline,"%lf,%lf,%lf,%lf,%lf",&r,&Q[i],&Qsca,&Qabs,&Gsca); // fills the Q[] array
			i++;
		}
		miefile.close();
		cout << "Processed mie file " << fname << "\n";
	}
	else cout << "\nUnable to open file " << fname << "\n";
	return;
}
void nadjust(double n[], int len, int offset, double f)
// adjusts the number density, n(r), with solution vector f
{
	int i;
	for(i=offset;i<len+offset;i++){
		n[i]*=f;
	}
	return;
}
void initn(double r[], double n[])
// initializes the number density, n(r).  The initial distribution function used should not affect the final solution.
{
	int i;
	for(i=0;i<619;i++){
		n[i]=N0/log(10.0)/r[i]*0.001*pow(r[i]*CONV_FACT,-4.0);  // one of Grassl's theoretical aerosol distributions
	}
	return;
}
void initH(matrix &H)
// initializes the smoothing matrix
{
	int i,j;
	for(i=0;i<6;i++){
		if(i==0){
			H[0][0]=1.0;
			H[0][1]=-2.0;
			H[0][2]=1.0;
			H[8][6]=1.0;
			H[8][7]=-2.0;
			H[8][8]=1.0;
		}
		else if(i==1){
			H[1][0]=-2.0;
			H[1][1]=5.0;
			H[1][2]=-4.0;
			H[1][3]=1.0;
			H[7][5]=1.0;
			H[7][6]=-4.0;
			H[7][7]=5.0;
			H[7][8]=-2.0;
		}
		else{
			H[i][i-2]=1.0;
			H[i][i-1]=-4.0;
			H[i][i]=6.0;
			H[i][i+1]=-4.0;
			H[i][i+2]=1.0;
		}
	}
	return;
}
void mtrxprnt(matrix &a)
// used for troubleshooting to print out a matrix 
{
	int i,j,row,col;
	row=a.size();
	col=a[0].size();
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			cout << a[i][j] << "\t";
		}
		cout<<endl;
	}
	return;
}
void mtrxmult(matrix &X,matrix &Y,matrix &Z)
// multiplies matrices: Z = X * Y
// X, Y, and Z are defined as matrices prior to call
{
	int i,j,k,Xrow,Xcol,Ycol;
	Xrow=X.size();
	Xcol=X[0].size();
	Ycol=Y[0].size();

//	cout<<"Initialize Z[i][j]..."<<endl;
	for(i=0;i<Xrow;i++){
		for(j=0;j<Ycol;j++){
			Z[i][j]=0.0;
		}
	}
//	mtrxprnt(Z);
//	cout<<endl<<"Calculate..."<<endl;
	for(i=0;i<Xrow;i++){
		for(j=0;j<Ycol;j++){
			for(k=0;k<Xcol;k++){
//				cout<<Z[i][j]<<"\t";
				Z[i][j]+=X[i][k]*Y[k][j];
			}
//			cout<<endl;
		}
	}
	return;
}

void Tmtrxmult(matrix &X,matrix &Y,matrix &Z)
// multiplies the transpose of a matrix with a matrix: Z = X^T * Y
// X, Y, and Z are defined as matrices prior to call
{
	int i,j,k,Xrow,Xcol,Ycol;
	Xrow=X.size();
	Xcol=X[0].size();
	Ycol=Y[0].size();

	for(i=0;i<Xcol;i++){
		for(j=0;j<Ycol;j++){
			Z[i][j]=0.0;
		}
	}
	for(i=0;i<Xcol;i++){
		for(j=0;j<Ycol;j++){
			for(k=0;k<Xrow;k++){
				Z[i][j]+=X[k][i]*Y[k][j];
			}
		}
	}
	return;
}

void Sinvmtrx(matrix &a, matrix &b)
/* Simplified invert matrix routine assumes non-zero diagonal elements
Input matrix a[1..N][1..N] is destroyed. Output (inverted) matrix is b[1..N][1..N]
*/
{
	int i,j,k,rows=a.size(),cols=a[0].size();
	double temp;

	// first make matrix b the identity matrix
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++){
			if(i==j){b[i][j]=1.0;}
			else{b[i][j]=0.0;}
		}
	}

	// begin the forward elimination
	for(i=0;i<rows;i++){
		temp=a[i][i];
		for(j=0;j<cols;j++){		// normalize diagonal element
			a[i][j] /= temp;	// and adjust all values to right
			b[i][j] /= temp;
		}
		if(i<rows-1){
			for(k=i+1;k<rows;k++){
				temp=a[k][i];
				for(j=0;j<cols;j++){
					a[k][j] -= temp*a[i][j];
					b[k][j] -= temp*b[i][j];
				}
			}
		}
	}
	// now do the back substitution
	for(i=rows-1;i>=0;i--){
		for(k=i-1;k>=0;k--){
			temp=a[k][i];
			for(j=0;j<cols;j++){
				a[k][j] -= temp*a[i][j];
				b[k][j] -= temp*b[i][j];
			}
		}
	}
	return;
}

double optimum_gamma(matrix &AtCiA, matrix &H, matrix &AtCi, double aod[])
// Optimizes the Lagrangian multiplier gamma for the matrix inversion routine
{
	matrix B(11,vector<double>(11,0.0));
	matrix Bi(11,vector<double>(11,0.0));
	matrix BiAtCi(11,vector<double>(6,0.0));
	double f[11],gamma,gammarel;
	int i,j,k;
	bool flag=false;
	
	// Loop to find minimum gammarel for which all f are positive
	k=1;
	while(k<=1000 && flag==false){
		// This section calculates the solution vector f
		gammarel=k/1000.0;
		gamma=gammarel*AtCiA[0][0]/H[0][0];
		for(i=0;i<11;i++){for(j=0;j<11;j++){B[i][j]=AtCiA[i][j]+gamma*H[i][j];}} // Calculate B matrix
		Sinvmtrx(B,Bi); // Calculate inverse B matrix
		mtrxmult(Bi,AtCi,BiAtCi);  // Calculate B inverse * A transposed * C inverse
		for(i=0;i<11;i++){f[i]=0.0;} // reset solution vector f
		for(i=0;i<11;i++){for(j=0;j<6;j++){f[i]+=BiAtCi[i][j]*aod[j];}} // Calculate solution vector f
		// test for negative (non-physical) solutions
		if(f[0]<0 || f[1]<0 || f[2]<0 || f[3]<0 || f[4]<0 || f[5]<0 || f[6]<0 || f[7]<0 || f[8]<0 || f[9]<0 || f[10]<0)  				{flag=false;}
		else {flag=true;}
		k++;
	}
//	cout << "k = " << k << endl;
//	for(i=0;i<11;i++){cout << f[i] << "\t";}
//	cout << endl << "gamma = " << gamma << endl;
	return (gamma);
}
