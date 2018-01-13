//===================================================
//Using two motif files to calculate motif-motif similarity
//Designed by Shaoqiang Zhang (zhangshaoqiang@gmail.com)
//SPIC.cpp
//====================================================
#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
using namespace std;
//-----------------------------------------------------
typedef vector<vector<double> > Matrix;
double MtfPrfSim(const Matrix& prfA, const Matrix& prfB);
//-----------------------------------------------------
int main(int argc, const char** argv){
    double cutoff;
    if((argc!=3)&&(argc!=4)){
        cout<<"\nCalculate motif-motif Similarity between Positions with Information Contents (SPIC).\n\nUSAGE:\n\n";
        cout<<"PROGRAM_NAME <Motif_A_file> <Motif_B_file> [OPTION: cutoff (-1<x<1)]\n\n";
        cout<<"NOTE: The format of a motif file must be as follows:\n";
        cout<<"(The lines with labels 'A', 'C', 'G', and 'T' form the profile matrix)\n";
        cout<<"(The line with label 'I' forms the vector of coloumn information contents)\n";
        cout<<"(The lines with labels 'a', 'c', 'g', and 't' form the frequency matrix)\n\n";
        cout<<"Motif_label\n";
        cout<<"A xx xx xx xx xx ... xx\n";
        cout<<"C xx xx xx xx xx ... xx\n";
        cout<<"G xx xx xx xx xx ... xx\n";
        cout<<"T xx xx xx xx xx ... xx\n";
        cout<<"I xx xx xx xx xx ... xx\n";
        cout<<"a xx xx xx xx xx ... xx\n";
        cout<<"c xx xx xx xx xx ... xx\n";
        cout<<"g xx xx xx xx xx ... xx\n";
        cout<<"t xx xx xx xx xx ... xx\n\n";
	cout<<"If you want to use the metric, please cite the paper:\n";
	cout<<"S. Zhang et al. A Novel Information Contents Based Similarity Metric for Comparing TFBS Motifs, Proceedings of ISB2012: 32-36.\n\n";
	cout<<"Designed by Shaoqiang Zhang (zhangshaoqiang@gmail.com), July 2012\n\n";
        exit(1);
    }
    ifstream motifAfile(argv[1]);
    ifstream motifBfile(argv[2]);
    if(argc==3){
        cutoff=-1;
    }else{
        if(argc==4){
            string scutoff(argv[3]);
            istringstream iscutoff(scutoff);
            iscutoff>>cutoff;
        }else{
            cout<<"Error input arguments, please read USAGE by only typing program name\n\n";
            exit(1);
        }
    }
    //--------------------------------------------------
    string Amotiflabel;
    string strA; getline(motifAfile,strA);
    istringstream AfirstLine(strA);
    AfirstLine>>Amotiflabel;
    Matrix Aprofile;
    char ACGT;
    for(string s;getline(motifAfile,s);){
        vector<double> Aline;
        istringstream sin(s);
        sin>>ACGT;
        for(double a; sin>>a; ){
            Aline.push_back(a);
        }
        Aprofile.push_back(Aline);
        Aline.clear();
    }
    string Bmotiflabel;
    string strB; getline(motifBfile,strB);
    istringstream BfirstLine(strB);
    BfirstLine>>Bmotiflabel;
    Matrix Bprofile;
    for(string t;getline(motifBfile,t);){
        vector<double> Bline;
        istringstream tin(t);
        tin>>ACGT;
        for(double b;tin>>b; ){
            Bline.push_back(b);
        }
        Bprofile.push_back(Bline);
        Bline.clear();
    }
    //-------------------------------------------
    double ScoreAB=MtfPrfSim(Aprofile,Bprofile);
    double ScoreBA=MtfPrfSim(Bprofile,Aprofile);
    double ScoreAA=MtfPrfSim(Aprofile,Aprofile);
    double ScoreBB=MtfPrfSim(Bprofile,Bprofile);
    double simScore1;
    (ScoreAB>ScoreBA)?simScore1=ScoreAB :simScore1=ScoreBA;
    double simScore2;
    (ScoreAA>ScoreBB)?simScore2=ScoreAA :simScore2=ScoreBB;
    double simScore3=simScore1/simScore2;
    double simValue;
    (simScore3>1)?simValue=1 :simValue=simScore3;
    if(simValue>=cutoff){
    cout.precision(2);
        cout<<Amotiflabel<<" "<<Bmotiflabel<<" "<<simValue<<endl;
    }
}
//---------------------------------------------------------------------------------------------------
double MtfPrfSim(const Matrix& prfA, const Matrix& prfB){
    //using the frequence matrix of motif A to scan the profile matrix of motif B
    int maxPositive=0;
    int countPositive=0;
    //double maxUnderScore=0.0;
    double tempscore=0.0;
    double maxScore=-1000.0;
    int mtflenA=prfA[0].size();
    int mtflenB=prfB[0].size();
    double seqnumA=prfA[5][0]+prfA[6][0]+prfA[7][0]+prfA[8][0];
    //cout<<mtflenA<<" "<<mtflenB<<" "<<seqnumA<<endl;
    //for(int i=21;i<22;++i){
    for(int i=0;i<mtflenA+mtflenB-1;++i){
        int lowB,uppB,lowA,uppA;
        (i-mtflenA+1>0)? lowB=i-mtflenA+1 : lowB=0;
        (i<mtflenB-1)? uppB=i : uppB=mtflenB-1;
        (mtflenA-i-1>0)? lowA=mtflenA-i-1 : lowA=0;
        (mtflenA-1<mtflenA+mtflenB-i-2)? uppA=mtflenA-1 : uppA=mtflenA+mtflenB-i-2;
	//cout<<lowB<<" "<<uppB<<" "<<lowA<<" "<<uppA<<endl;
        Matrix frequenceAmatch, profileBmatch;
        for (int ib=0; ib<5; ++ib){
            vector<double> Bvector;
            for(int bb=lowB;bb<=uppB;++bb){
                Bvector.push_back(prfB[ib][bb]);
            }
            profileBmatch.push_back(Bvector);
            //Bvector.clear();
        }
        for (int ia=5;ia<9;++ia){
            vector<double> Avector;
            for(int aa=lowA;aa<=uppA;++aa){
                Avector.push_back(prfA[ia][aa]);
            }
            frequenceAmatch.push_back(Avector);
            //Avector.clear();
        }
        countPositive=0;
        tempscore=0.0;
        for (int nn=0;nn<=uppB-lowB;++nn){
            double compareColumn=(profileBmatch[0][nn])*(frequenceAmatch[0][nn]);
            for(int mm=1;mm<4;++mm){
                compareColumn=compareColumn+(profileBmatch[mm][nn])*(frequenceAmatch[mm][nn]);
            }
            double columnProfile=(profileBmatch[4][nn])*compareColumn;
            tempscore=tempscore+columnProfile;
            if(columnProfile>=0){
                countPositive++;
            }
        }
        if(countPositive>=maxPositive){
            maxPositive=countPositive;
        //}
            if(tempscore>maxScore){
                maxScore=tempscore;
            }
        }
    }
    return maxScore;
}
