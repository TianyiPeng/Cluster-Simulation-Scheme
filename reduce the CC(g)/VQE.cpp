#include<iostream>
#include<cstring>
#include<cstdlib>
#include<cstdio>
#include<ctime>
#include<complex>
#include<cassert>
#include "Eigen/Dense"

#define pow2(n) (1<<(n))
#define rep(i,r) for (int i = 0; i < r; i++)
using namespace std;
using namespace Eigen;


const int n = 6;//number of qubits
const int M = 1<<n;//Hilbert space dimension
const int Small_Depth = 5;//sign for Small Depth algorithm 
const int Small_Depth_Value = 3;//the depth value of Small Depth algorithm
const int Depth = 9;//depth for the general algorithms
const int tot_comp = 6;//the number of different algorithms

const int hamil_num = 50;//the number of Pauli term 

int low_enable[tot_comp][Depth+1];//control which entanglers are enable
int ham[hamil_num][n];//store the Pauli term of the Hamiltonian
double weight[hamil_num];//store the weight of each Pauli term

double ave[tot_comp][300008];//record the result for averaging
int total_test = 100000;//number of iterations
double num_exp = 100;//number of experiments

struct quantum_gate//1-qubit gate class
{
    complex<double> gate[2][2];
    friend quantum_gate operator*(const quantum_gate &a,
            const quantum_gate &b)//matrix multiplication
            {
                quantum_gate G;
                rep(i,2)
                    rep(j,2)
                    {
                        G.gate[i][j] = 0;
                        rep(k,2)
                            G.gate[i][j] += a.gate[i][k]*b.gate[k][j];
                    }
                return G;
            }
}Pauli[4];
complex<double> tmp[M];
struct quantum_state
{
    complex<double> a[M];
    void initial()//initial the state to |00..0>
    {
        for (int i = 0; i < M; i++)
            a[i] = 0;
        a[0] = 1;    
    }
    void CNOT(int st, int ed)//Apply CNOT gate 
    //control qubit: st, target qubit: ed
    {
        int l1 = min(st, ed);
        int l2 = max(st, ed);
        assert(l1 != l2);
        for (int s1=0; s1<pow2(l1); s1++)
            for (int s2=0;s2<pow2(l2-l1-1); s2++)
                for (int s3=0; s3<pow2(n-l2-1); s3++)
                //enumerate the state of the qubits other than st and ed
                {  
                    int s = (s3 << (l2+1)) + (s2 << (l1+1)) + s1 + (1<<st);
                    swap(a[s], a[s+(1<<ed)]);
                }
    }
    void gate(quantum_gate U, int ed)//apply single-qubit gate to ed
    {
        for (int s1 = 0; s1 < pow2(ed); s1++)
            for (int s2=0; s2 < pow2(n-ed-1); s2++)
            {
                int s = (s2 << (ed+1)) + s1;
                complex<double> x1 = a[s];
                complex<double> x2 = a[s+(1<<ed)];
                a[s] = U.gate[0][0] * x1 + U.gate[0][1] * x2;
                a[s+(1<<ed)] = U.gate[1][0] * x1 + U.gate[1][1] * x2;
            }
    }
    double estimation(int ham_term[])//
    {
        complex<double> ans = 0;
        rep(i,M)
            tmp[i] = a[i];
        rep(i,n)
        {
            gate(Pauli[ham_term[i]], i);
        }
        rep(i,M)
        {
            ans += a[i]*conj(tmp[i]);
            a[i] = tmp[i];
        }
        assert(abs(ans.imag())<1e-8);
        return ans.real();
    }
    double estimation()//
    {
        double ans = 0;
        rep(i,hamil_num)
            ans += estimation(ham[i])*weight[i];
        return ans;
    }
};
int seed;
void initial_problem()
{
    seed = time(0);
    srand(seed);
    rep(i,hamil_num)
        rep(j,n)
        {
            ham[i][j] = (rand()&3);
            weight[i] = ((double)rand())/RAND_MAX*2-1;
        }
}

quantum_gate X_gate(double beta)//generate the exp(-i*beta*X) gate
{
    quantum_gate U;
    complex<double> i(0.0,1.0); 
    U.gate[0][0] = (exp(-i*beta) + exp(i*beta))/2.0;
    U.gate[1][1] = (exp(-i*beta) + exp(i*beta))/2.0;
    U.gate[0][1] = (exp(-i*beta) - exp(i*beta))/2.0;
    U.gate[1][0] = (exp(-i*beta) - exp(i*beta))/2.0;
    return U;
} 
quantum_gate Z_gate(double beta)//generate the exp(-i*beta*Z) gate
{
    quantum_gate U;
    complex<double> i(0.0,1.0); 
    U.gate[0][0] = exp(-i*beta);
    U.gate[1][1] = exp(i*beta);
    U.gate[0][1] = 0;
    U.gate[1][0] = 0;
    return U;
} 
quantum_gate ZXZ_gate(double b1, double b2, double b3)
{
    return Z_gate(b3)*X_gate(b2)*Z_gate(b1);   
}

quantum_state Q;
double estimate(double beta[], int S)//build the circuit
{
        Q.initial();
        int t = 0;
        rep(i,n)
        {
            Q.gate(ZXZ_gate(0, beta[t], beta[t+1]), i);
            t += 2; 
        }
        int D = Depth;
        if (S==Small_Depth)
            D = Small_Depth_Value;
        for (int i=1; i<= D; i++)
        {
            for (int j = 0; j < n-1; j++)
                for (int k = j+1; k < n; k++)
                if ((j+1)*2 <= n && (k+1)*2 > n)   
                {
                    if (low_enable[S][i])
                        Q.CNOT(j,k);
                }
                else
                    Q.CNOT(j, k);
            for (int j = 0; j < n; j++)
            {
                Q.gate(ZXZ_gate(beta[t], beta[t+1], beta[t+2]), j);
                t += 3;
            }
        }
        return Q.estimation();
}


complex<double> H[M][M];
double find_classical_value()//find the optimal value of classically
{
    MatrixXcf A(M,M);
    rep(sx,M)
        rep(sy,M)
            A(sx,sy) = 0;
    rep(i,hamil_num)
    {
        rep(sx,M)
            rep(sy,M)
            {
                complex<double> tmp = weight[i];
                rep(k,n)
                {
                    tmp *= Pauli[ham[i][k]].gate[(sx>>k)&1][(sy>>k)&1];
                }
                A(sx,sy) += tmp;
            }
    }
    SelfAdjointEigenSolver<MatrixXcf> eigensolver(A);
    if (eigensolver.info() != Success) abort();
    cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << endl;
    return eigensolver.eigenvalues().minCoeff();
}

void quantum_optimization_SPSA()
{
    int L = (n)*(3*Depth+2);
    double beta[tot_comp][L], d[L], delta[L];
    //initial_value
    for (int i = 0; i < L; i++)
        rep(j,tot_comp)
            beta[j][i] = 1;
    int T = 0;
    double an;
    double cn;
    double F[tot_comp];

    double vc = find_classical_value();
    rep(test,total_test)
    {
        T++;
        //set parameters
        an = 0.3 / pow(T, 0.6);
        cn = 0.1 / pow(T, 0.1);
        //SPSA
        rep(j,tot_comp)
        {
            for (int i = 0; i < L; i++)
                delta[i] = (rand()&1)*2-1;
            for (int i = 0; i < L; i++)
                beta[j][i] += delta[i]*cn;
            double F_plus = estimate(beta[j], j);
            for (int i = 0; i < L; i++)
                beta[j][i] -= delta[i]*cn*2;
            double F_minus = estimate(beta[j], j);
            for (int i = 0; i < L; i++)
                beta[j][i] += delta[i]*cn;
            for (int i = 0; i < L; i++)
            {
                beta[j][i] -= an * (F_plus-F_minus)/(2*cn*delta[i]);
            }
        }
        rep(j,tot_comp)
        {
            F[j] =  estimate(beta[j], j);
            ave[j][test] += abs((F[j]-vc)/vc);
        }
        if (T % 1000 == 0)
        {
            printf("%d %d ", T, seed);
            rep(k,tot_comp)
                printf("%lf ", abs((F[k]-vc)/vc));
            printf("\n");
        }
    }     
    cout << vc << endl;
}

int main()
{
    Pauli[0].gate[0][0] = Pauli[0].gate[1][1] = 1;// I gate
    Pauli[1].gate[0][1] = Pauli[1].gate[1][0] = 1;// X gate
    Pauli[2].gate[0][0] = 1; Pauli[2].gate[1][1] = -1;// Z gate
    Pauli[3].gate[0][1] = -complex<double>(0, 1); 
    Pauli[3].gate[1][0] = -Pauli[3].gate[0][1];//Y gate

    rep(i,Depth) low_enable[0][i+1] = low_enable[Small_Depth][i+1] = 1;
    low_enable[1][3] = low_enable[1][5] = low_enable[1][7] = 1;
    low_enable[2][3] = low_enable[2][7] = 1;
    low_enable[3][5] = 1;
    
    rep(i,num_exp)
    {
        initial_problem();
        quantum_optimization_SPSA();
        cout << i << endl;
    }
    
    freopen("result.txt", "w", stdout);
    rep(i,total_test)
    {
        printf("%d ", i);
        rep(j,tot_comp)
            printf("%lf ", ave[j][i]/num_exp);
        printf("\n");
    }
    
    return 0;
}