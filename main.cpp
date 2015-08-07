#include<iostream>
#include<algorithm>
#include<stdlib.h>
#include<cstddef>
#include<stdio.h>

using namespace std;

int maximum(int x,int y, int z)
{
    return max({x,y,z});
}
int main()
{
    char yes,yn;
    do
    {
        string sequence1,sequence2;
        string protein_sequence="ARNDCQEGHILKMFPSTWYV";
        size_t n1,n2;
        size_t mm1,mm2;
        int gap= -2, gap_extend= -8;

        int BLOSUM50[][20]=
        {
            {5,-2,-1,-2,-1,-1,-1,0,-2,-1,-2,-1,-1,-3,-1,1,0,-3,-2,0},
            {-2,7,-1,-2,-4,1,0,-3,0,-4,-3,3,-2,-3,-3,-1,-1,-3,-1,-3},
            {-1,-1,7,2,-2,0,0,0,1,-3,-4,0,-2,-4,-2,1,0,-4,-2,-3},
            {-2,-2,2,8,-4,0,2,-1,-1,-4,-4,-1,-4,-5,-1,0,-1,-5,-3,-4},
            {-1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1},
            {-1,1,0,0,-3,7,2,-2,1,-3,-2,2,0,-4,-1,0,-1,-1,-1,-3},
            {-1,0,0,2,-3,2,6,-3,0,-4,-3,1,-2,-3,-1,-1,-1,-3,-2,-3},
            {0,-3,0,-1,-3,-2,-3,8,-2,-4,-4,-2,-3,-4,-2,0,-2,-3,-3,-4},
            {-2,0,1,-1,-3,1,0,-2,10,-4,-3,0,-1,-1,-2,-1,-2,-3,2,-4},
            {-1,-4,-3,-4,-2,-3,-4,-4,-4,5,2,-3,2,0,-3,-3,-1,-3,-1,4},
            {-2,-3,-4,-4,-2,-2,-3,-4,-3,2,5,-3,3,1,-4,-3,-1,-2,-1,1},
            {-1,3,0,-1,-3,2,1,-2,0,-3,-3,6,-2,-4,-1,0,-1,-3,-2,-3},
            {-1,-2,-2,-4,-2,0,-2,-3,-1,2,3,-2,7,0,-3,-2,-1,-1,0,1},
            {-3,-3,-4,-5,-2,-4,-3,-4,-1,0,1,-4,0,8,-4,-3,-2,1,4,-1},
            {-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3},
            {1,-1,1,0,-1,0,-1,0,-1,-3,-3,0,-2,-3,-1,5,2,-4,-2,-2},
            {0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,2,5,-3,-2,0},
            {-3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1,1,-4,-4,-4,15,2,-3},
            {-2,-1,-2,-3,-3,-1,-2,-3,2,-1,-1,-2,0,4,-3,-2,-2,2,8,-1},
            {0,-3,-3,-4,-1,-3,-3,-4,-4,4,1,-3,1,-1,-3,-2,0,-3,-1,5}
        };
        system("cls");
        bool invalid;
        do
        {
            invalid=false;
            cout<<"Enter Sequence1 ( only from these alphabets: "<<protein_sequence<<" ): ";
            cin>>sequence1;
            cout<<"Enter Sequence2 ( only from these alphabets: "<<protein_sequence<<" ): ";
            cin>>sequence2;
            for (string::iterator i = sequence1.begin(); sequence1.end() != i; i++)
            {
                *i = toupper(*i);
            }
            for (string::iterator i = sequence2.begin(); sequence2.end() != i; i++)
            {
                *i = toupper(*i);
            }
            mm1 = sequence1.find_first_not_of(protein_sequence);
            mm2 = sequence2.find_first_not_of(protein_sequence);
            if(mm1!=string::npos)
            {
                cout<<endl<<"Invalid characters in Protein sequence: "<<sequence1[mm1]<<" at the position: "<<mm1+1<<endl;
                invalid=true;
            }
            else if(mm2!=string::npos)
            {
                cout<<endl<<"Invalid characters in Protein sequence"<<sequence2[mm2]<<" at the position: "<<mm2+1<<endl;
                invalid=true;
            }
        }
        while(invalid);
        system("cls");
        int length1=sequence1.length();
        int length2=sequence2.length();
        cout<<endl<<"Sequences entered by you are:"<<endl<<endl;
        cout<<"Sequence 1:\t"<<sequence1<<endl
            <<"Sequence 2:\t"<<sequence2<<endl;
        int S[length2+1][length1+1];
        int i,j;
        int MX[length2+1][length1+1];
        int MY[length2+1][length1+1];
        int  M[length2+1][length1+1];
        int pointer[length2+1][length1+1];


        /**
        ********************************************
         Initializing all matrix.
        ********************************************
        **/

        for(i=0; i<length2+1; i++)
        {
            for(j=0; j<length1+1; j++)
            {
                M[i][j] = 0;
                MX[i][j] = 0;
                MY[i][j] = 0;
                S[i][j] = 0;
                pointer[i][j] = 0;
            }
        }
        for(j=0; j<length1+1; j++)
        {
            if(j==0)
            {
                continue;
            }
            else if(j==1)
            {
                S[0][j]=j*gap;
            }
            else
            {
                S[0][j]= S[0][j-1] + gap_extend;
            }
        }
        for(i=0; i<length2+1; i++)
        {
            if(i==0)
            {
                continue;
            }
            else if(i==1)
            {
                S[i][0]=i*gap;
            }
            else
            {
                S[i][0]=S[i-1][0] + gap_extend;
            }
        }


        /**
        ********************************************
          Parent Matrix M[i][j], MX[i][j], MY[i][j]
        ********************************************
        */

        for(i=1; i<length2+1; i++)
        {
            for(j=1; j<length1+1; j++)
            {
                int IXM = M[i-1][j] + gap;
                int IIX = MX[i-1][j] + gap_extend;
                MX[i][j] = max({IXM, IIX});

                int IYM = M[i][j-1] + gap;
                int IIY = MY[i][j-1] + gap_extend;
                MY[i][j] = max({IYM, IIY});

                n1=protein_sequence.find(sequence1[j-1]);
                n2=protein_sequence.find(sequence2[i-1]);
                int MM = M[i-1][j-1] + BLOSUM50[n1][n2];
                int MIX = MX[i-1][j-1] + BLOSUM50[n1][n2];
                int MIY = MY[i-1][j-1] + BLOSUM50[n1][n2];
                M[i][j] = max({MM, MIX, MIY});
                S[i][j] = max({M[i][j], MX[i][j], MY[i][j]});


                if(S[i][j] == 0)
                {
                    pointer[i][j] = 0; //0 means end of the path
                    continue;
                }
                if(S[i][j] == MX[i][j])
                {
                    pointer[i][j] = 1; //1 means trace up
                    continue;

                }
                else if(S[i][j] == MY[i][j])
                {
                    pointer[i][j] = 2; //2 means trace left
                    continue;
                }
                else if(S[i][j] == M[i][j])
                {
                    pointer[i][j] = 3; // 3 means trace diagonal
                    continue;
                }
            }
        }

        /**
        *****************************************************
         Printing all three parent matrix and pointer matrix
        *****************************************************
        **/
/*
        cout<<endl<<"Matrix M"<<endl<<endl;
        for(i=0; i<length2+1; i++)
        {
            for(j=0; j<length1+1; j++)
            {
                printf("%5d ",M[i][j]);
            }
            cout<<endl;
        }
        cout<<endl<<"Matrix MX"<<endl<<endl;
        for(i=0; i<length2+1; i++)
        {
            for(j=0; j<length1+1; j++)
            {
                printf("%5d ",MX[i][j]);
            }
            cout<<endl;
        }
        cout<<endl<<"Matrix MY"<<endl<<endl;
        for(i=0; i<length2+1; i++)
        {
            for(j=0; j<length1+1; j++)
            {
                printf("%5d ",MY[i][j]);
            }
            cout<<endl;
        }
        cout<<endl<<"Pointer Matrix"<<endl<<endl;
        for(i=0; i<length2+1; i++)
        {
            for(j=0; j<length1+1; j++)
            {
                printf("%5d ",pointer[i][j]);
            }
            cout<<endl;
        }

*/

        /**
        ******************************************
          displaying scoring matrix
        ******************************************
        */


        cout<<endl<<"The Dynamic Programming Matrix:"<<endl<<endl;
        cout<<"\t";
        for(unsigned int i=0; i<=sequence1.size(); i++)
        {
            if(i==0)
            {
                cout<<"   ";
            }
            else
            {
                cout<<"    "<<sequence1[i-1]<<" ";
            }
        }
        cout<<endl<<endl;
        for(i=0; i<length2+1; i++)
        {
            if(i==0)
            {
                cout<<"     ";
            }
            else
            {
                cout<<"   "<<sequence2[i-1]<<" ";
            }
            for(int j=0; j<length1+1; j++)
            {
                printf("%5d ",S[i][j]);
            }
            cout<<endl<<endl;
        }


        /**
        ******************************************
          trace back process
        ******************************************
        */

        i=length2;
        j=length1;
        int score=S[i][j];
        string aln1="";
        string aln2="";

        while(i>0 && j>0)
        {
            if(pointer[i][j]==3)
            {
                aln1+=sequence1[j-1];
                aln2+=sequence2[i-1];
                i--;
                j--;
            }
            else if(pointer[i][j]==2)
            {
                aln1+=sequence1[j-1];
                aln2+="_";
                j--;
            }
            else if(pointer[i][j]==1)
            {
                aln1+="_";
                aln2+=sequence2[i-1];
                i--;
            }
        }
        if(j>0)
        {
            while(j>0)
            {
                aln1+=sequence1[j-1];
                aln2+="_";
                j--;
            }
        }
        else if(i>0)
        {
            while(i>0)
            {
                aln1+="_";
                aln2+=sequence2[i-1];
                i--;
            }
        }
        reverse(aln1.begin(),aln1.end());
        reverse(aln2.begin(),aln2.end());
        cout<<endl<<endl<<"Aligned Sequences are:"<<endl<<endl<<aln1<<endl<<aln2<<endl;
        cout<<endl<<"Optimum alignment score: "<<score<<endl<<endl;

        /**
        ******************************************
          Exit
        ******************************************
        */


        do
        {
            cout<<"Do again? (Y/N):";
            cin>>yes;
            yn = toupper(yes);
            if(yn!='N' && yn != 'Y')
            {
                cout<<endl<<"Invalid selection"<<endl;
            }
            if(yn=='N' || yn == 'Y')
            {
                break;
            }

        }
        while(true);
    }
    while(yn == 'Y');
}

