// DigitRecog_RealTime.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <conio.h>
using namespace std;
int counter = 1, accuracy = 0;
long double pi_avg[10][6], a_avg[10][6][6], b_avg[10][6][33];
void checkfile(fstream* f, string fname){
	if(!f) cout<<"Unable to open file "<<fname<<endl;
	//else   cout<<"Successfully opened file "<<fname<<endl;
}
long double GetDCShift()
{						
	int x,i=0;
	long double DCShift = 0;
	long double max=0;
	fstream fin("Noise.txt");	checkfile(&fin,"Noise.txt");			//Open and read file "Noise.txt"
	while(fin>>x){
		DCShift += x;	
		i++;
		}
	DCShift /= i;
	return DCShift; 
}
long double DCShift = GetDCShift();
long double codebook[33][13];
void readcodebook(){
	fstream fin("codebook.txt"); checkfile(&fin,"codebook.txt");
	for(int i = 1; i<=32; i++){
		for(int j = 1; j<=12; j++){
			fin>>codebook[i][j];
		}
	}
}

class computeCi{
	#define N 320						//put framesize to be considered
	#define p 12
	#define pi 3.1415926535
public: 
	long double s[321],R[p+1],a[p+1],c[p+1],RSW_c[p+1]; 
	long double Normalize;
	computeCi(int frame[],long double ci[]){
		Normalize = NormFactor(frame);
		process_frame(frame);
		HammingWindow(frame);
		Calc_Ri(frame);
		Calc_ai();
		Calc_ci();
		Calc_ci_RSW();
		for(int i = 1; i<=p; i++){
			ci[i] = RSW_c[i];
		}
	}
long double NormFactor(int frame[])
{
	long double maxx = 0,minn = 0,limit = 5000, scale;
	for(int i = 1; i<=320; i++){ 
			maxx = (frame[i]>maxx)? frame[i] : maxx;
			minn = (frame[i]<minn)? frame[i] : minn;			
		}
	scale = (maxx - minn) / 2;
	return (limit/scale);
}
void process_frame(int frame[]){
	for(int i = 1; i<=320; i++){
		frame[i] = (frame[i]-DCShift) * Normalize;
	}
}
void HammingWindow(int frame[])      //cap the frame with hamming window
{
	long double w[N];
	for(int m=0; m<=N-1; m++){
		w[m] = 0.54 - 0.46 * cos((double)(2*22*m)/(7*(N-1))) ;
	}
	for(int i = 1; i<=N; i++){
		frame[i] = frame[i] * w[i-1];	
	}
}
void Calc_Ri(int frame[])		///returns Ri s of the given frame
{
	for(int k = 0; k<=p; k++)
	{
		long double temp = 0;
		for(int y = 1; y<= N-k; y++)
		{
			temp += frame[y]*frame[y+k];
		}
		R[k] = temp;
	}
	
}
void Calc_ai(){  // impl Levinsion Alg
	long double E[p+1],K[p+1];
	long double A[p+1][p+1];
	E[0] = R[0];
	for(int i=1; i<=p; i++){
		long double temp = 0;
		for(int j = 1; j<= i-1; j++){
			temp += A[j][i-1] * R[i-j];
		}
		K[i] = (R[i] - temp) / E[i-1];
		A[i][i] = K[i];
		for(int j = 1; j<= i-1; j++){
			A[j][i] = A[j][i-1] - K[i] * A[i-j][i-1];
		}
		E[i] = (1 - K[i]*K[i]) * E[i-1];
	}
	for(int i = 1; i<=p; i++){
		a[i] = A[i][p];
	}
}
void Calc_ci(){  //compute cepstral coeff
	long double energy = R[0];
	c[0] = log(energy*energy);
	for(int m = 1; m<=p; m++)
	{
		long double temp = 0;
		for(int k = 1; k<=m-1;k++)
			temp+= ((double)k/m) * c[k] * a[m-k];
		c[m] = a[m] +temp;
	}	
}
void Calc_ci_RSW(){  //compute cepstral coeff with raised sine window
	long double energy = R[0],Cweight[p];
	RSW_c[0] = log(energy*energy);
	for(int m = 1; m<=p; m++)
	{
		long double temp = 0;
		for(int k = 1; k<=m-1;k++)
			temp+= ((double)k/m) * RSW_c[k] * a[m-k];
		RSW_c[m] = a[m] + temp;
	}
	for(int m = 0; m<p; m++)
		Cweight[m] = 1 + (p/2)*sin((double)(22*m)/(7*p));

	for(int i = 1; i<=p; i++)
		RSW_c[i] *= Cweight[i-1];
}
};
long double tokhura(long double c[], long double cr[]){  //returns tokhura distance
	long double d=0;
	long double TW[13] = {0,1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0};
	for(int i = 1; i<=12; i++){
		d+= TW[i] * ((c[i] - cr[i])*(c[i] - cr[i]));
	}
	return d;
	}
int findvq(long double ci[]){
	int vq = 0;
	long double tokhuradistance, minn = LDBL_MAX;
	for(int i = 1; i<=32; i++){
		tokhuradistance = tokhura(ci,codebook[i]);
		vq = (tokhuradistance < minn) ? i : vq;
		minn = (tokhuradistance < minn) ? tokhuradistance : minn;
	}
	return vq;
}
void initialmodel(long double Pi[6], long double a[6][6], long double b[6][33], int d){
	if(counter == 1){
				fstream fpi("pi.txt"); checkfile(&fpi,"pi.txt");
				for(int i = 1; i<=5; i++){
					fpi>>Pi[i];
				}
				fpi.close();
				fstream fa("a.txt"); checkfile(&fa,"a.txt");
				for(int i = 1; i<=5; i++){
					for(int j = 1; j<=5; j++){
						fa>>a[i][j];
					}
				}
				fa.close();
				fstream fb("b.txt"); checkfile(&fb,"b.txt");
				for(int i = 1; i<=5; i++){
					for(int j = 1; j<=32; j++){
						fb>>b[i][j];
					}
				}
				fb.close();
	}else{
				string pifile = "./Final Model/" + to_string((long long) d) + "_" + "pi.txt"; 
				string afile = "./Final Model/" + to_string((long long) d) + "_" + "a.txt"; 
				string bfile = "./Final Model/" + to_string((long long) d) + "_" + "b.txt"; 
				fstream fpi(pifile); checkfile(&fpi,pifile);
				for(int i = 1; i<=5; i++){
					fpi>>Pi[i];
				}
				fpi.close();

				fstream fa(afile); checkfile(&fa,afile);
				for(int i = 1; i<=5; i++){
					for(int j = 1; j<=5; j++){
						fa>>a[i][j];
					}
				}
				fa.close();
				fstream fb(bfile); checkfile(&fb,bfile);
				for(int i = 1; i<=5; i++){
					for(int j = 1; j<=32; j++){
						fb>>b[i][j];
					}
				}
				fb.close();	
	}
	
}
class createmodel{
public:
	long double a[6][6], b[6][33], pie[6], alpha[151][6], beta[151][6], probability,
				delta[151][6],p_star, last_val_of_pStar,
				a_bar[6][6], b_bar[6][33], pie_bar[6], zai[151][6][6], gamma[151][6];
	int T,O[150],psi[151][6],q_star[151];
	createmodel(int o[], int T, long double matpi[], long double mata[6][6], long double matb[6][33], fstream &hout, int digit){
		initialization(o,T,matpi,mata,matb);
		int itr = 1;
		do{
		last_val_of_pStar = p_star;
		validate(); manipulateBmat(); validate();
		forward(); backward(); viterbi(); 
		reestimate();update();
		}while(p_star > last_val_of_pStar && (itr < 100 || p_star - last_val_of_pStar > 1e-50)); 
		manipulateBmat(); validate();
		printModel(hout);
		findavg(digit);
		}
	void initialization(int o[], int T, long double matpi[], long double mata[6][6], long double matb[6][33]){
		for(int i = 1; i<=T; i++){
			O[i] = o[i];
		}
		for(int i = 1; i<=5; i++){
			for(int j = 1; j<=5; j++){
				a[i][j] = mata[i][j];
			}
		}
		for(int i = 1; i<=5; i++){
			for(int j = 1; j<=32; j++){
				b[i][j] = matb[i][j];
			}
		}
		for(int i = 1; i<=5; i++){
			pie[i] = matpi[i];
		}
		this->T = T;
		for(int t = 0; t<=150; t++){
			for(int i = 0; i<=5; i++){
				gamma[t][i] = alpha[t][i] = beta[t][i] = delta[t][i] = 0;
			}
		}

	}
	void forward(){
		for(int i = 1; i<=5; i++){
			int x = O[1];
			alpha[1][i] = pie[i] * b[i][O[1]];
		}
		for(int t = 1; t<=T-1; t++){
			for(int j = 1; j<=5; j++){
					long double temp = 0;
					for(int i = 1; i<=5; i++){
						temp += alpha[t][i] * a[i][j];
					}
					int x = O[t+1];
					alpha[t+1][j] = temp * b[j][O[t+1]];
					//printf("%Lf ",alpha[t+1][j]);//<<" ";
			}
			//cout<<endl;
		}
		probability = 0;
		for(int i = 1; i<=5; i++){
			probability += alpha[T][i];
		}
	}
	void backward(){
		for(int i=1; i<=5; i++){
			beta[T][i] = 1;
		}
		for(int t=T-1; t>=1; t--){
			for(int i=1; i<=5; i++){
				long double temp = 0;
				for(int j = 1; j<=5; j++){
					temp += a[i][j] * b[j][O[t+1]] * beta[t+1][j];
				}
				beta[t][i] = temp;
			}
		}
	}
	void viterbi(){
		for(int i = 1; i<=5; i++){
			delta[1][i] = pie[i] * b[i][O[1]];
			psi[1][i] = 0;
		}
		for(int t = 2; t<=T; t++){
			for(int j = 1; j<=5; j++){
				int maxi = 1;
				long double maxx = delta[t-1][1] * a[1][j];
				for(int i = 1; i<=5; i++){
					long double temp = delta[t-1][i] * a[i][j];
					maxi = temp > maxx ? i : maxi;
					maxx = temp > maxx ? temp : maxx;
				}
				delta[t][j] = maxx * b[j][O[t]];
				psi[t][j] = maxi;
			}
		}
		p_star = delta[T][1];
		q_star[T] = 1;
		for(int i = 1; i<=5; i++){
			q_star[T] = delta[T][i] > p_star ? i : q_star[T];
			p_star = delta[T][i] > p_star ? delta[T][i] : p_star;
		}
		for(int t = T-1; t>=1; t--){
			q_star[t] = psi[t+1][q_star[t+1]];
		}
	}
	void reestimate(){
		for(int t = 1; t<=T-1; t++){
			for(int i = 1; i<=5; i++){
				for(int j = 1; j<=5; j++){
					long double temp = alpha[t][i] * a[i][j] * b[j][O[t+1]] * beta[t+1][j];
					zai[t][i][j] = temp / probability;
					//printf("%Lf ",zai[t][i][j]);
				}
			}
			//cout<<endl;
		}
		//cout<<endl;
		for(int t = 1; t<=T-1; t++){
			for(int i = 1; i<=5; i++){
				long double temp = 0;
				for(int j = 1; j<=5; j++){
					temp += zai[t][i][j];
				}
				gamma[t][i] = temp;
				//printf("%Lf ",gamma[t][i]);
			}
			//cout<<endl;
		}
		//cout<<endl;
		for(int i =1; i<=5; i++){
			pie_bar[i] = gamma[1][i];
			//printf("%Lf ",pie_bar[i]);
		}
		//cout<<endl<<endl;
		for(int i =1; i<=5; i++){
			for(int j =1; j<=5; j++){
				long double numr = 0 , denr = 0;
				for(int t = 1; t<=T-1; t++){
					numr += zai[t][i][j];
					denr += gamma[t][i];
				}
				a_bar[i][j] = numr / denr;
				//printf("%Lf ",a_bar[i][j]);
			}
			//cout<<endl;
		}
		//cout<<endl;
		for(int j = 1; j<=5; j++){
			for(int k = 1; k<=32; k++){
				long double denr = 0.0, numr = 0.0;
				for(int t = 1; t<=T; t++){
					numr += (O[t] == k) ? gamma[t][j] : 0.0;
				}

				for(int t = 1; t<=T; t++){
					denr += gamma[t][j];
					//printf("%Lf ",denr);
				}

				b_bar[j][k] = ((long double) numr) / denr;
				//printf("%Lf ",denr);
				//cout<<endl;
			}
			//cout<<endl;
		}
	}
	void update(){
		for(int i = 1; i<=5; i++){
			pie[i] = pie_bar[i];
			for(int j = 1; j<=5;j++){
				a[i][j] = a_bar[i][j];
			}
			for(int k = 1; k<=32; k++){
				b[i][k] = b_bar[i][k];
				//printf("%Lf ",b[i][k]);//<<" ";
			}
			//cout<<endl;
		}
	}
	void printModel(fstream &hout){
		for(int i = 1; i<=5; i++){
			hout<<pie[i]<<" ";
		}
		hout<<endl<<endl;
		for(int i = 1; i<=5; i++){
			for(int j = 1; j<=5; j++){
				hout<<a[i][j]<<" ";
			}
			hout<<endl;
		}
		hout<<endl<<endl;
		for(int i = 1; i<=5; i++){
			for(int j = 1; j<=32; j++){
				hout<<b[i][j]<<" ";
			}
			hout<<endl;
		}
		hout<<endl<<endl;
		hout<<p_star<<endl<<endl;
		for(int i=1; i<=T; i++){
			hout<<O[i]<<" ";
		}
		hout<<endl;
		for(int i=1; i<=T; i++){
			hout<<q_star[i]<<" ";
		}
	}
	void validate(){
		long double sum = 0;
		for(int i = 1; i<=5; i++){
			sum += pie[i];
		}
	
		for(int i = 1; i<=5; i++){
			sum = 0; int maxi = 1;
			for(int j = 1; j<=5; j++){
				sum += a[i][j];
				maxi = (a[i][j] > a[i][maxi]) ? j : maxi;
			}
			a[i][maxi] += (1-sum>0)? 1 - sum : 0;
			a[i][maxi] -= (sum-1>0)? sum - 1 : 0;
		}

		for(int i = 1; i<=5; i++){
			sum = 0; int maxi = 1;
			for(int j = 1; j<=32; j++){
				sum += b[i][j];
				maxi = (b[i][j] > a[i][maxi]) ? j : maxi;
			}
			b[i][maxi] += (1-sum>0)? 1 - sum : 0;
			b[i][maxi] -= (sum-1>0)? sum - 1 : 0;
		}

	}
	void manipulateBmat(){
		for(int i = 1; i<=5; i++){
			for(int j = 1; j<=32; j++){
				b[i][j] = (b[i][j] == 0) ? 1e-40 : b[i][j];
			}
		}
	}
	void findavg(int d){
			for(int i = 1; i<=5; i++){			
				pi_avg[d][i] += pie[i];
				for(int j = 1; j<=5; j++){
					a_avg[d][i][j] += a[i][j];
				}
				for(int k = 1; k<=32; k++){
					b_avg[d][i][k] += b[i][k];
				}
			}
		
	}
};
void trainmodel(){
	fstream ufile; ufile.open("Universe.csv",fstream::out); checkfile(&ufile,"Universe.csv");
	readcodebook();
	
	while(counter<2){
		long double Pi[6], a[6][6], b[6][33];
		for(int d = 0; d<=9; d++){
			for(int i = 0; i<=5; i++){			/*initialize pi_avg, a_avg, b_avg to 0*/
				pi_avg[d][i] = 0;
				for(int j = 0; j<=5; j++){
					a_avg[d][i][j] = 0;
				}
				for(int k = 0; k<=32; k++){
					b_avg[d][i][k] = 0;
				}
			}
		}
	for(int d = 0; d<=9; d++){
		initialmodel(Pi,a,b,1);
		for(int u = 1; u<=20; u++){
			string recordedFile = "./Recorded Files/234101011_E_" + to_string((long long) d) + "_" + to_string((long long) u) + ".txt";
			string cepstralFile = "./Cepstrals/234101011_E_" + to_string((long long) d) + "_" + to_string((long long) u) + ".txt";
			string vqFile = "./Observation Sequence/234101011_E_" + to_string((long long) d) + "_" + to_string((long long) u) + ".txt";
			string logfile = "./HMM logs/234101011_E_" + to_string((long long) d) + "_" + to_string((long long) u) + ".txt";

			fstream fin(recordedFile); checkfile(&fin,recordedFile);							//read recorded files
			int speech[100000] , x , len = 0;
			while(fin>>x){
				len++;
				speech[len] = x;
			}
			fin.close();

			int frame[151][321]; long double ci[151][13]; int fno = 1;			//compute ci
			for(int i = 0; i + 320 <= len && fno <= 150; i+=100, fno++){
				for(int j = 1; j <= 320; j++){
					frame[fno][j] = speech[i+j];
				}
			}
			int no_of_frames = fno - 1;  
			for(int i = 1; i<=no_of_frames; i++){
				computeCi Ciobj = computeCi(frame[i],ci[i]);
			}
			fstream fout; fout.open(cepstralFile,fstream::out); checkfile(&fout,cepstralFile);
			for(int i = 1; i<=no_of_frames; i++){
				for(int j = 1; j<=12; j++){
					//cout<<ci[i][j]<<" ";
					fout<<ci[i][j]<<" ";
					ufile<<ci[i][j]<<",";
				}
				//cout<<endl;
				fout<<endl;
				ufile<<endl;
			}
			fout.close();

			int vq[151];
			fstream vout; vout.open(vqFile,fstream::out); checkfile(&vout,vqFile); 
			for(int i = 1; i<=no_of_frames; i++){
				vq[i] = findvq(ci[i]);
				vout<<vq[i]<<" ";
			}
			vout.close();
			
			fstream hout; hout.open(logfile,fstream::out); checkfile(&hout,logfile);
			createmodel hmm = createmodel(vq,no_of_frames,Pi,a,b,hout,d);
			hout.close();
		}
	}
		for(int d = 0; d<=9; d++){
				for(int i = 0; i<=5; i++){			/*initialize pi_avg, a_avg, b_avg to 0*/
					pi_avg[d][i] /= 20;
					for(int j = 0; j<=5; j++){
						a_avg[d][i][j] /= 20;
					}
					for(int k = 0; k<=32; k++){
						b_avg[d][i][k] /= 20;
					}
				}
			}
	for(int d = 0; d<=9; d++){
			string pifile = "./Final Model/" + to_string((long long) d) + "_" + "pi.txt"; 
			string afile = "./Final Model/" + to_string((long long) d) + "_" + "a.txt"; 
			string bfile = "./Final Model/" + to_string((long long) d) + "_" + "b.txt"; 
			fstream piout; piout.open(pifile,fstream::out); checkfile(&piout,pifile);
			fstream aout; aout.open(afile,fstream::out); checkfile(&aout,afile);
			fstream bout; bout.open(bfile,fstream::out); checkfile(&bout,bfile);

			for(int i = 1; i<=5; i++){			
				piout<<pi_avg[d][i]<<" "; 
				for(int j = 1; j<=5; j++){
					aout<<a_avg[d][i][j]<<" ";
				}
				for(int k = 1; k<=32; k++){
					bout<<b_avg[d][i][k]<<" ";
				}
				aout<<endl; bout<<endl;
			}
	}
	counter++;
	}

}
long double call_forward(int O[],int T,long double pie[6], long double a[6][6], long double b[6][33]){
	long double alpha[151][6];
	for(int t = 0; t<=150; t++){
		for(int i=0; i<=5; i++){
			alpha[t][i]=0;
		}
	}
	for(int i = 1; i<=5; i++){
			int x = O[1];
			alpha[1][i] = pie[i] * b[i][O[1]];
		}
		for(int t = 1; t<=T-1; t++){
			for(int j = 1; j<=5; j++){
					long double temp = 0;
					for(int i = 1; i<=5; i++){
						temp += alpha[t][i] * a[i][j];
					}
					int x = O[t+1];
					alpha[t+1][j] = temp * b[j][O[t+1]];
					//printf("%Lf ",alpha[t+1][j]);//<<" ";
			}
			//cout<<endl;
		}
		long double probability = 0;
		for(int i = 1; i<=5; i++){
			probability += alpha[T][i];
		}
		return probability;
	}
void testmodel(){
	readcodebook();
	for(int d = 0; d<=9; d++){
			string pifile = "./Final Model/" + to_string((long long) d) + "_" + "pi.txt"; 
			string afile = "./Final Model/" + to_string((long long) d) + "_" + "a.txt"; 
			string bfile = "./Final Model/" + to_string((long long) d) + "_" + "b.txt"; 
			fstream piin; piin.open(pifile); checkfile(&piin,pifile);
			fstream ain; ain.open(afile); checkfile(&ain,afile);
			fstream bin; bin.open(bfile); checkfile(&bin,bfile);

			for(int i = 1; i<=5; i++){			
				piin>>pi_avg[d][i]; 
				for(int j = 1; j<=5; j++){
					ain>>a_avg[d][i][j];
				}
				for(int k = 1; k<=32; k++){
					bin>>b_avg[d][i][k];
				}
			}
	}

	system("Recording_Module.exe 4 input_file.wav input_file.txt"); 

	cout<<"press any key\n"; 
			string recordedFile = "input_file.txt";
			string cepstralFile = "Cepstral_Value.txt";
			string vqFile = "Observation_Sequence.txt";
			
			fstream fin(recordedFile); checkfile(&fin,recordedFile);							//read recorded files
			int speech[100000] , x , len = 0;
			while(fin>>x){
				len++;
				speech[len] = x;
			}
			fin.close();

			int frame[151][321]; long double ci[151][13]; int fno = 1;			//compute ci
			for(int i = 0; i + 320 <= len && fno <= 150; i+=100, fno++){
				for(int j = 1; j <= 320; j++){
					frame[fno][j] = speech[i+j];
				}
			}
			int no_of_frames = fno - 1;  
			for(int i = 1; i<=no_of_frames; i++){
				computeCi Ciobj = computeCi(frame[i],ci[i]);
			}
			fstream fout; fout.open(cepstralFile,fstream::out); checkfile(&fout,cepstralFile);
			for(int i = 1; i<=no_of_frames; i++){
				for(int j = 1; j<=12; j++){
					//cout<<ci[i][j]<<" ";
					fout<<ci[i][j]<<" ";
				}
				//cout<<endl;
				fout<<endl;
			}
			fout.close();

			int vq[151];
			fstream vout; vout.open(vqFile,fstream::out); checkfile(&vout,vqFile); 
			for(int i = 1; i<=no_of_frames; i++){
				vq[i] = findvq(ci[i]);
				cout<<vq[i]<<" ";
				vout<<vq[i]<<" ";
			}
			vout.close();
			cout<<"\nTesting for file"<<recordedFile<<endl;
			int  maxi = 0;
			long double maxx = -1, probability;
			for(int digit = 0; digit <=9; digit++){
				probability = call_forward(vq,no_of_frames,pi_avg[digit],a_avg[digit],b_avg[digit]);
				cout<<"Probability(O / lambda"<<digit<<") = "<<probability<<endl;
				maxi = (probability>maxx) ? digit : maxi;
				maxx = (probability>maxx) ? probability : maxx;
			}
			cout<<"Digit recognized by HMM : "<<maxi<<endl;
			cout<<"--------------------------------------------------------------------------------------------------------------------------------------------------------\n";
		
}
int _tmain(int argc, _TCHAR* argv[])
{	
	
	//trainmodel();	cout<<"Training Done"<<endl;		/* Dont run Training and testing at once.*/
	testmodel();

	getch();
	return 0;
}

