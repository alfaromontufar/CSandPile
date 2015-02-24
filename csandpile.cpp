#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <set>
#include <vector>
#include <list>

using namespace std;

ofstream fout;

class csandpile{
	
private:
	
	int *degree;
	int *configuration;
	int *config;
	int **rest;
	int **laplacian;
	double **inverse;
	int m;
	int sink;
	list< vector<int> > configs;
	
public:
	
	csandpile(char *);
	int wini(int *);
	void sum(int, int *);
	void top(int *);
	void top(int *, int*);
	void print();
	void print(int *);
	void print(float *);
	void print(double *);
	void printmatrix(int **, int);
	void printmatrix(float **,int);
	void stable();
	bool areequals(int *, int *);
	void recurrent(int *);
	int powerof(int *);
	void powerof();
	void identity();
	void reprec();
	void reprecinv();
	int cofactor(int **, int);
	void determinant();
	int inversematrix();
	void minimaldegree(int *);
	void minimaldegree();
	void group();
	void all();
	void table();
};

csandpile::csandpile(char *name){
	
	
	char *input_file_name = new char [(strlen(name)+3)];
	strcpy(input_file_name,name);
	strcat(input_file_name,".gph");  
	ifstream fin(input_file_name);
	
	int dim;
	
	fin>>dim;
	
	laplacian = new int*[dim];
	for(int i=0; i<dim; i++)
		laplacian[i] = new int [dim];  
	
	for(int i=0; i<dim; i++)
		for(int j=0; j<dim; j++)
			fin>>laplacian[i][j]; 
	
	fin>>sink;
	
	m = dim - 1;
	
	rest = new int*[m];
	for(int i=0; i<m; i++)
		rest[i] = new int [m];  
	
	int flagi=0, flagj=0;
	for(int i=0; i<dim; i++)
		for(int j=0; j<dim; j++){
			if(i!=(sink-1))
				if(j!=(sink-1))
					rest[i-flagi][j-flagj]=laplacian[i][j];
				else
					flagj=1;
				else
					flagi=1;
			if(j==(dim-1))
				flagj=0;
		}
	
	degree = new int [m];
	
	for(int i=0; i<m; i++)
		degree[i] = rest[i][i]; 
	
	config = new int[dim];
	configuration = new int[m];
	
	int count = 0;
	
	while(!fin.eof()){
		vector<int> v(m);
		for(int i=0; i<dim; i++)
			fin>>config[i];
		flagi=0;
		for(int i=0; i<dim; i++)
			if(i!=(sink-1)){
				if(count == 0)
					configuration[i-flagi] = config[i];
				v[i-flagi] = config[i];
			}
			else
				flagi=1;
		configs.push_front(v);
		if(count == 0)
			count = 1;
	}
	
	if(count == 0)
		for(int i=0; i<m; i++)
			configuration[i]=degree[i]-1;
	
	fin.close();
}

int csandpile::wini(int *tocheck){  // Where is not it?
	for(int i=0; i<m; i++)
		if(tocheck[i]>=degree[i])
			return i;
	return -1;
}

void csandpile::sum(int k, int *tocheck){
	for(int i=0; i<m; i++)
		tocheck[i] -=  rest[k][i];
}

void csandpile::top(int *tocheck){
	int aux=wini(tocheck);
	
	while(aux!=-1){
		sum(aux, tocheck);
		aux=wini(tocheck);
	}
}

void csandpile::top(int *tocheck, int *times){
	for(int i=0; i<m; i++)
		times[i] = 0;
	
	int aux=wini(tocheck);
	
	while(aux!=-1){
		times[aux]++;
		sum(aux, tocheck);
		aux=wini(tocheck);
	}
}

void csandpile::print(){                 //imprime configuration
	fout<<" ";
	for(int i=0; i<m; i++){
		if(i==(sink-1))
			fout<<"S ";
		fout<<configuration[i]<<" ";
	}
	if(m==(sink-1))
		fout<<"S";
	fout<<endl;
}

void csandpile::print(int *tocheck){                 //imprime configuration
	fout<<" ";
	for(int i=0; i<m; i++){
		if(i==(sink-1))
			fout<<"S ";
		fout<<tocheck[i]<<" ";
	}
	if(m==(sink-1))
		fout<<"S";
	fout<<endl;
}

void csandpile::print(float *tocheck){                 //imprime configuration
	fout<<" ";
	for(int i=0; i<m; i++){
		if(i==(sink-1))
			fout<<"S ";
		fout<<tocheck[i]<<" ";
	}
	if(m==(sink-1))
		fout<<"S";
	fout<<endl;
}

void csandpile::print(double *tocheck){                 //imprime configuration
	fout<<" ";
	for(int i=0; i<m; i++){
		if(i==(sink-1))
			fout<<"S ";
		fout<<tocheck[i]<<" ";
	}
	if(m==(sink-1))
		fout<<"S";
	fout<<endl;
}

void csandpile::printmatrix(int **matrix, int n){    //imprime matrix
	fout<<endl;
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)
			fout<<"  "<<matrix[i][j];
		fout<<endl;
	}
}

void csandpile::printmatrix(float **matrix, int n){  //imprime matrix
	fout<<endl;
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)
			fout<<"  "<<matrix[i][j];
		fout<<endl;
	}
}

void csandpile::stable(){
	int *aux;
	
	aux = new int [m];
    
	for(int i=0; i<m; i++)
		aux[i]=configuration[i];
	
	fout<<" The stable configuration of"<<endl;
	print(aux);
	top(aux);
	fout<<" is"<<endl;
	print(aux);
}

bool csandpile::areequals(int *aux, int * tocheck){
	for(int i=0; i<m; i++)
		if(aux[i]!=tocheck[i])
			return false;
	
	return true;
}

void csandpile::recurrent(int *tocheck){
	int *aux = new int [m];
	
	for(int i=0; i<m; i++)
		aux[i] = 2*(degree[i]-1);
	
	top(aux);
	
	for(int i=0; i<m; i++)
		aux[i] = 2*(degree[i]-1) - aux[i] + tocheck[i];
	
	top(aux);
	
	for(int i=0; i<m; i++)
		tocheck[i] = aux[i];
}

int csandpile::powerof(int * tocheck){
	int *aux1, *aux2;
	
	aux1 = new int [m];
	aux2 = new int [m];
	
	for(int i=0; i<m; i++)
		aux1[i] = tocheck[i];
	
	recurrent(aux1);
	
	for(int i=0; i<m; i++)
		aux2[i] = aux1[i];
	
	fout<<" Checking configuration: ";
	print(aux1);
	fout<<endl;
	
	fout<<" Powers"<<endl;
	
	int j=1;
	
	fout<<j++<<" - ";
	print(aux2);
	
	for(int i=0; i<m; i++)
		aux2[i] += aux1[i];
	
	int * times = new int[m];
	
	top(aux2,times);
	
	while(!areequals(aux1,aux2)){
		fout<<j++<<" - ";
		print(aux2);
		fout<<" T: ";
		print(times);
		for(int i=0; i<m; i++)
			aux2[i] += aux1[i];
		
		top(aux2,times);
	}  
	
	return (--j);
}

void csandpile::powerof(){
	powerof(configuration);
}

void csandpile::identity(){
	int *aux = new int [m];
	
	for(int i=0; i<m; i++)
		aux[i] = 2*(degree[i]-1);
	
	top(aux);
	
	for(int i=0; i<m; i++)
		aux[i] = 2*(degree[i]-1) - aux[i];
	
	top(aux);
	
	print(aux);
}

void csandpile::reprec(){
	int *aux = new int[m];
	
	for(int i=0; i<m; i++)
		aux[i] = configuration[i];
	
	recurrent(aux);
	
	print(aux);
}

void csandpile::reprecinv(){
	
	int *aux = new int [m];
	
	int max=1;
	for(int i=0; i<m; i++)
		if((degree[i]-1)!=0)
			if(max<((int) configuration[i]/(degree[i]-1)+1))
				max = (int) configuration[i]/(degree[i]-1)+1;
	
	max = max + 2;
	
	for(int i=0; i<m; i++)
		aux[i] = max*(degree[i]-1);
	
	top(aux);
	
	for(int i=0; i<m; i++)
		aux[i] = max*(degree[i]-1) - aux[i] - configuration[i];
	
	top(aux);
	
	print(aux);
}

int csandpile::cofactor(int **matrix, int n){
	
	if(n == 1)
		return matrix[0][0];
	else if(n == 2)
		return (matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1]);
	
	int sum=0;
	int **temp = new int *[n-1];
	for(int i=0; i<(n-1); i++)
		temp[i] = new int [n-1];
	
	for(int i=0; i<n; i++){
		for(int j=1; j<n; j++)
			for(int k=0; k<n; k++)
				if( k<i )
					temp[j-1][k] = matrix[j][k];
				else if( k>i )
					temp[j-1][k-1] = matrix[j][k];	 	  	      
		if(matrix[0][i]!=0)
			sum += matrix[0][i]* (int) pow(-1,i)*cofactor(temp,n-1);
	}
	return sum;
}

void csandpile::determinant(){  
	fout<<" Determinant: "<< cofactor(rest,m)<<endl;
}

int csandpile::inversematrix(){
	
	inverse = new double *[m];
	double **temp = new double *[m];
	
	for(int i=0; i<m; i++)
		inverse[i] = new double [m];
	
	for(int i=0; i<m; i++)
		temp[i] = new double [m];
	
	for(int i=0; i<m; i++)
		for(int j=0; j<m; j++)
			inverse[i][j] = 0;
	
	for(int i=0; i<m; i++)
		inverse[i][i] = 1;
	
	for(int i=0; i<m; i++)
		for(int j=0; j<m; j++)
			temp[i][j]=rest[i][j];
	
	for(int i=0; i<m; i++){
		int nzero=-1;
		for(int j=i; j<m; j++)
			if(temp[j][i]!=0){
				nzero=j;
				break;
			}
		if(nzero!=-1){
			if(nzero!=i){
				for(int j=0; j<m; j++){
					double aux=temp[i][j], iaux=inverse[i][j];
					temp[i][j] = temp[nzero][j];
					temp[nzero][j] = aux;
					inverse[i][j] = inverse[nzero][j];
					inverse[nzero][j] = iaux;
				}
			}
			
			double piv = temp[i][i];
			for(int j=0; j<m; j++){
				temp[i][j] = temp[i][j]/piv;
				inverse[i][j] = inverse[i][j]/piv;
			}
			
			for(int j=0; j<m; j++)
				if(i!=j && temp[j][i]!=0){
					double aux = temp[j][i];
					for(int k=0; k<m; k++){
						temp[j][k] = temp[j][k] - aux*temp[i][k];
						inverse[j][k] = inverse[j][k] - aux*inverse[i][k];
					}
				}
		}
		else
			return 0;
	}	
	return 1;
}

int gcd(int a, int b){
	
	int auxa=a, auxb=b, auxc;
	
	if(auxa<0)
		auxa=-auxa;
	
	if(auxb<0)
		auxb=-auxb;
	
	if (auxa<auxb){
		auxc=auxb;
		auxb=auxa;
		auxa=auxc;
	}
	
	while (auxb>0){
		auxc=auxa%auxb;
		auxa=auxb;
		auxb=auxc;
	} 
	
	return auxa;
	
}

void csandpile::minimaldegree(int * tocheck){
	double * min = new double [m];
	
	for(int i=0; i<m; i++)
		min[i]=0;
	
	if(inversematrix()){
		for(int i=0; i<m; i++)  
			for(int j=0; j<m; j++)  
				min[i] = min[i] + (double) inverse[i][j]*tocheck[j];
		
		for(int i=0; i<m; i++){
			min[i]-= (int) min[i];
			int aux=1;
			if(min[i] > 0.0)
				while((min[i] - (int) min[i]) != 0.0){
					aux = 10*aux;
					min[i] = 10*min[i];
				}
			min[i] = aux/gcd(aux,(int) min[i]);
		}
		
		int deg= (int) min[0];
		for(int i=1; i<m; i++)
			deg = gcd(deg, (int) min[i]);
		fout<<" The minimal degree of "<<endl;
		print(min);
		fout<<" is "<<endl<<" "<<deg;
	}
	else
		fout<<" Error: The reduced Laplacian Matrix is not invertible";
}

void csandpile::minimaldegree(){
	minimaldegree(configuration);
}

void csandpile::group(){
	
	int det = cofactor(rest,m);
	
	fout<<endl<<" Number of elements: "<<det<<endl;	
	
	int ngen = m;
	
	int ** generators = new int *[ngen];
	
	int * powers = new int[m];
	
	for(int i=0; i<ngen; i++)
		generators[i] = new int [m];
	
	for(int i=0; i<ngen; i++)
		for(int j=0; j<m; j++)
			generators[i][j] = 0;
	
	for(int i=0; i<ngen; i++){
		generators[i][i]=1;
		
		fout<<" Generator "<<(i+1)<<": ";       
		print(generators[i]);
		
		recurrent(generators[i]);
		
		powers[i] = powerof(generators[i]);
	}
	
	long mult=1;
	for(int i=0; i<ngen; i++)
		mult = mult * (long) powers[i];
	
	fout<<endl<<" Number of combinations: "<<mult;
	
	// it work with 1000000000
	if(mult<10000000){
		
		int * div = new int[ngen];
		
		for(int i=0; i<ngen; i++)
			div[i]=1;
		
		for(int i=0; i<ngen; i++)
			for(int j=0; j<ngen; j++)
				if(i > j)
					div[i] = div[i]*powers[j];
		
		//for(int i=0; i<nmin; i++)
		//  fout<<endl<<" div["<<i<<"]= "<<div[i];
		
		vector<int> v(m);
		set< vector<int> > conj;
		
		int *ident = new int [m];
		
		for(int i=0; i<m; i++)
			ident[i] = 2*(degree[i]-1);
		
		top(ident);
		
		for(int i=0; i<m; i++)
			ident[i] = 2*(degree[i]-1) - ident[i];
		
		top(ident);
		
		for(int j=0; j<m; j++)
			v[j]=ident[j];
		
		conj.insert(v);
		
		for(double i=1; i<mult; i++){
			int * sum = new int [m];
			for(int j=0; j<m; j++)
				sum[j]=0;
			
			for(int j=0; j<ngen; j++)
				for(int k=0; k<m; k++)
					sum[k]+=(  (int) (i/ (double) div[j]) % powers[j] )*generators[j][k];	
			
			top(sum);
			
			for(int j=0; j<m; j++)
				v[j]=sum[j];
			
			conj.insert(v);
			//print(sum);
		}
		
		fout<<endl<<" Size of all diferent linear combinations: "<<conj.size()<<endl<<endl;
		
	}
	else
		fout<<endl<<" It will not serve";
}

void csandpile::all(){
	
	int det = cofactor(rest,m);
	
	fout<<endl<<" Number of elements: "<<det<<endl;
	
	int *sequence = new int[m];
	int *adjvert = new int[m];
	int *value = new int[m];
	int nmin=0;    
	
	for(int i=0; i<(m+1); i++){ //first we know adjcent vertices of the sink
		if( i > (sink-1) ){
			if( laplacian[sink-1][i] < 0 )
				adjvert[i-1] = 1;
			else
				adjvert[i-1] = 0;
		}
		else if( i < (sink-1) )
			if( laplacian[sink-1][i] < 0 )
				adjvert[i] = 1;
			else
				adjvert[i] = 0;
	}    
	
	for(int i=0; i<m; i++)
		nmin+=adjvert[i];
	
	for(int i=0; i<m; i++){ //we give the first value
		value[i]=degree[i];
		if(adjvert[i] == 1)
			value[i]--;
	}
	
	int **generators = new int *[nmin];
	int *powers = new int [nmin];
	
	for(int i=0; i<nmin; i++)
		generators[i] = new int [m];
	
	int auxi=0;
	for(int i=0; i<m; i++){
		if(adjvert[i] == 1){
			
			int *marked = new int[m];
			int *auxvalue = new int[m];
			
			for(int j=0; j<m; j++){
				marked[j] = 0;
				auxvalue[j] = value[j];
			}
			
			sequence[i]=0;
			marked[i] = 1;
			
			for(int j=0; j<m; j++)
				if(i != j)
					auxvalue[j]+=rest[i][j];
			
			
			for(int j=1; j<m; j++){
				int findex;
				for(int  k=0; k<m; k++)	 
					if(marked[k] == 0){
						if(adjvert[k] == 1)
							findex = k;
						else
							for(int l=0; l<m; l++)
								if(marked[l] == 1 && rest[k][l] < 0){
									findex = k;
									break;
								}
						break;
					}
				
				sequence[findex] = j;
				marked[findex] = 1;
				for(int k=0; k<m; k++)
					if(marked[k] == 0)
						auxvalue[k]+=rest[findex][k];
			}
			
			fout<<endl<<" Sequence: "<<endl;
			for(int j=0; j<m; j++){
				generators[auxi][j] = auxvalue[j];
				fout<<" "<<(sequence[j]+1);
			}
			fout<<endl;
			powers[auxi] = powerof(auxvalue);
			auxi++;
		}
	}
	
	//Here we compute the linear combinations of all generators
	
	long mult=1;
	for(int i=0; i<nmin; i++)
		mult = mult * (long) powers[i];
	
	fout<<endl<<" Number of combinations: "<<mult;
	
	if(mult<1000000000){
		
		int * div = new int[nmin];
		
		for(int i=0; i<nmin; i++)
			div[i]=1;
		
		for(int i=0; i<nmin; i++)
			for(int j=0; j<nmin; j++)
				if(i > j)
					div[i] = div[i]*powers[j];
		
		//for(int i=0; i<nmin; i++)
		//  fout<<endl<<" div["<<i<<"]= "<<div[i];
		
		vector<int> v(m);
		set< vector<int> > conj;
		
		int *ident = new int [m];
		
		for(int i=0; i<m; i++)
			ident[i] = 2*(degree[i]-1);
		
		top(ident);
		
		for(int i=0; i<m; i++)
			ident[i] = 2*(degree[i]-1) - ident[i];
		
		top(ident);
		
		for(int j=0; j<m; j++)
			v[j]=ident[j];
		
		conj.insert(v);
		
		for(double i=1; i<mult; i++){
			int * sum = new int [m];
			for(int j=0; j<m; j++)
				sum[j]=0;
			
			for(int j=0; j<nmin; j++)
				for(int k=0; k<m; k++)
					sum[k]+=(  (int) (i/ (double) div[j]) % powers[j] )*generators[j][k];	
			
			top(sum);
			
			for(int j=0; j<m; j++)
				v[j]=sum[j];
			
			conj.insert(v);
			//print(sum);
		}
		
		fout<<endl<<" Size of all diferent linear combinations: "<<conj.size()<<endl<<endl;
		
	}
	else
		fout<<endl<<" It will not serve";
	
}

void csandpile::table(){	
	
	int ngen = configs.size();
	
	int * powers = new int[m];
	
	list < list < vector<int> > > gens;
	
	for(int i=0; i<ngen; i++){
		
		list < vector<int> > auxl;
		
		vector<int> auxv(m);
		
		auxv = configs.back();
		
		configs.pop_back();
		
		int *aux = new int [m];
		
		for(int j=0; j<m; j++)
			aux[j] = auxv[j];
		
		fout<<" Generator "<<(i+1)<<": ";       
		print(aux);
		
		recurrent(aux);
		
		for(int k=0; k<m; k++)
			auxv[k] = aux[k];
		
		auxl.push_back(auxv);		
		
		int *aux1;
		
		aux1 = new int [m];
		
		for(int j=0; j<m; j++)
			aux1[j] = aux[j];
		
		fout<<" Checking configuration: ";
		print(aux1);
		fout<<endl;
		
		fout<<" Powers"<<endl;
		
		int j=1;
		
		fout<<j++<<" - ";
		print(aux);
		
		for(int k=0; k<m; k++)
			aux[k] += aux1[k];
		
		int * times = new int[m];
		
		top(aux,times);
		
		while(!areequals(aux1,aux)){
			
			fout<<j++<<" - ";
			print(aux);
			//fout<<" T: ";
			//print(times);
			
			for(int k=0; k<m; k++){
				auxv[k] = aux[k];
				aux[k] += aux1[k];
			}				
			
			auxl.push_back(auxv);
			
			top(aux,times);
		}  
		
		powers[i]=--j;
		
		gens.push_back(auxl);
	}
	
	long mult=1;
	for(int i=0; i<ngen; i++)
		mult = mult * (long) powers[i];
	
	fout<<endl<<" Number of combinations: "<<mult;
	
	// it work with 1000000000
	if(mult<10000){
		
		int * div = new int[ngen];
		
		for(int i=0; i<ngen; i++)
			div[i]=1;
		
		for(int i=0; i<ngen; i++)
			for(int j=0; j<ngen; j++)
				if(i > j)
					div[i] = div[i]*powers[j];
		
		vector<int> v(m);
		set< vector<int> > conj;
		
		conj.insert(gens.back().back()); //insert the identity
		
		fout<<endl<<" Linear combinations";
		
		for(double i=0; i<mult; i++){
			int * sum = new int [m];
			for(int j=0; j<m; j++)
				sum[j]=0;
			
			list< list<vector<int> > >::iterator itgens;
			list< vector<int> >::iterator itlist;
			
			itgens = gens.begin();
			
			fout<<endl;
			for(int j=0; j<ngen; j++){
				
				list<vector<int> > auxl;
				int indcons = ((int) (i/ (double) div[j]) % powers[j]);
				
				auxl = *itgens;
				itlist = auxl.begin();
				
				int l=0;
				
				while(itlist != auxl.end()){
					vector<int> auxv;
					auxv = *itlist;
					if(l == indcons){
						for(int k=0; k<m; k++)
							sum[k]+=auxv[k];
						break;
					}
					itlist++;
					l++;
				}
				
				fout<<" "<<(indcons+1);
				itgens++;
				
			}
			
			//fout<<" bfr ";
			//print(sum);
			
			int *times = new int [m];
			top(sum,times);
			
			for(int j=0; j<m; j++)
				v[j]=sum[j];
			
			conj.insert(v);
			//fout<<" Toppled ";
			fout<<"  - ";
			print(sum);
			//fout<<"          ";
			//print(times);
		}
		
		fout<<endl<<" Size of all diferent linear combinations: "<<conj.size()<<endl<<endl;
		
	}
	else
		fout<<endl<<" It will not serve";	
	
}


int main(int argc, char * argv[]){
	
	if(argc == 1){
		cout<<"No input file"<<endl;
	}
	
	else if(argc == 2){
		csandpile sand(argv[1]);
		char *output_file_name = new char [(strlen(argv[1])+3)];
		strcpy(output_file_name,argv[1]);
		strcat(output_file_name,".csp");
		fout.open(output_file_name, ios::trunc);
		fout<<endl<<" Write, after of the executable file, one of the following options:"<<endl;
		fout<<endl<<" -s           to obtain the stable configuration";
		fout<<endl<<" -p           to obtain the powers of the recurrent configuration";
		fout<<endl<<" -i           to obtain the identity";
		fout<<endl<<" -r           to obtain the recurrent configuration";
		fout<<endl<<" -ri          to obtain the inverse recurrent configuration";
		fout<<endl<<" -det         to obtain the determinant of the reduced Laplacian matrix";
		fout<<endl<<" -group       to obtain the powers of the standard base";
		fout<<endl<<" -sum         to obtain the sandpile sum of the two vectors";
		fout<<endl<<" -complete n  to create the Laplacian matrix of the complete graph of n vertices";
		fout<<endl<<" -path n      to create the Laplacian matrix of the path of n vertices";
		fout<<endl<<" -cycle n     to create the Laplacian matrix of the cycle of n vertices";
		fout<<endl<<" -version     to see the CSandPile version";
		fout<<endl<<endl;
		fout.close();
	}
	
	else if(argc==3){
		csandpile sand(argv[1]);
		char *output_file_name = new char [(strlen(argv[1])+3)];
		strcpy(output_file_name,argv[1]);
		strcat(output_file_name,".csp");
		fout.open(output_file_name, ios::trunc);
		
		if( strcmp( argv[2] , "-s") == 0 )
			sand.stable();
		else if( strcmp(argv[2] , "-p") == 0 )
			sand.powerof();
		else if( strcmp(argv[2] , "-i") == 0 ){
			fout<<" Identity: ";
			sand.identity();
		}
		else if( strcmp(argv[2] , "-r") == 0 ){
			fout<<" The recurrent configuration of ";
			sand.print();
			fout<<" is";      
			sand.reprec();
		}
		else if( strcmp(argv[2] , "-ri") == 0 ){
			fout<<" Inverse recurrent configuration: ";
			sand.reprecinv();
		}
		else if( strcmp(argv[2] , "-det") == 0 )
			sand.determinant();
		else if( strcmp(argv[2] , "-deg") == 0 )
			sand.minimaldegree();
		else if( strcmp(argv[2] , "-group") == 0 )
			sand.group();
		else if ( strcmp(argv[2] , "-all") == 0 )
			sand.all();
		else if ( strcmp(argv[2] , "-table") == 0 )
			sand.table();
		else if( strcmp(argv[2] , "-version") == 0 )
			fout<<" CSandPile version 0.5.10.6.24";
		else
			fout<<"Error: Not recognized command"<<endl;
		fout.close();
	}
	
	else if (argc==4){
		ofstream foutg;
		char *input_file_name = new char [(strlen(argv[1])+3)];
		strcpy(input_file_name,argv[1]);
		strcat(input_file_name,".gph");
		foutg.open(input_file_name, ios::trunc);
		int dim = atoi(argv[3]);
		
		if( strcmp(argv[2] , "-complete") == 0){      
			foutg<<dim<<endl;
			if(dim>0){
				for(int i=0; i<dim; i++){
					for(int j=0; j<dim; j++)
						if(i==j)
							foutg<<(dim-1)<<" ";
						else
							foutg<<"-1 ";
					foutg<<endl;
				}
				foutg<<dim;
			}
		}
		else if( strcmp(argv[2] , "-path") == 0){
			foutg<<dim<<endl;
			if(dim>0){
				for(int i=0; i<dim; i++){
					for(int j=0; j<dim; j++)
						if(i==j)
							if(i==0 || i==dim-1)
								foutg<<"1 ";
							else
								foutg<<"2 ";
							else
								if(i==(j+1) || i==(j-1))
									foutg<<"-1 ";
								else
									foutg<<"0 ";
					foutg<<endl;
				}
				foutg<<dim;
			}
		}
		else if( strcmp(argv[2] , "-cycle") == 0){
			foutg<<dim<<endl;
			if(dim>0){
				for(int i=0; i<dim; i++){
					for(int j=0; j<dim; j++)
						if(i==j)
							foutg<<"2 ";
						else if(i==(j+1) || i==(j-1))
							foutg<<"-1 ";
						else if((i==(dim-1))&&(j==0))
							foutg<<"-1 ";
						else if((j==(dim-1))&&(i==0))
							foutg<<"-1 ";
						else
							foutg<<"0 ";
					foutg<<endl;
				}
				foutg<<dim;
			}
		}
		else
			fout<<"Error: Not recognized command"<<endl;
		foutg.close();
	}  
	
	return 0;
}
