#include<RcppArmadillo.h>

using namespace Rcpp;

RcppExport SEXP placsolver1cpp(SEXP X, SEXP S, SEXP THETA, SEXP qrbeval, SEXP QRB, SEXP SAPV, SEXP NM, SEXP randoms1, SEXP randoms2) {

Rcpp::NumericMatrix Xr(X);                 
Rcpp::NumericVector Sr(S);
Rcpp::NumericMatrix THETA2r(THETA);
Rcpp::NumericMatrix qrbevalr(qrbeval);
Rcpp::NumericVector QRBr(QRB);
Rcpp::NumericVector SAPVr(SAPV);
Rcpp::NumericVector NMr(NM);
Rcpp::NumericMatrix randoms1r(randoms1);  
Rcpp::NumericMatrix randoms2r(randoms2);

int N0m = SAPVr.size();
int N0 = N0m + 1;
int B = THETA2r.nrow();
int N = NMr(0);
int M = NMr(1);
int Nm = N - 1;
int n = N*M;
int antony = THETA2r.ncol();
int blob = M*N0m;

arma::mat x(Xr.begin(), M, 2, false);
arma::vec s(Sr.begin(), N0, false);
arma::mat t(THETA2r.begin(), B, antony, false);
arma::mat QRBEVAL(qrbevalr.begin(), N, N0, false);
arma::vec qrb(QRBr.begin(), QRBr.size(), false);
arma::vec sapv(SAPVr.begin(), N0m, false);
arma::mat LUCY1(randoms1r.begin(), blob, B, false);
arma::mat LUCY2(randoms2r.begin(), blob, B, false);

double A1 = 7.5;
double A2 = 0;
double B1;
//double B2;

arma::vec ssapv(N0m);
ssapv = sqrt(sapv);

arma::uvec FF(N0m);
arma::uvec SS(N0m);
FF(0) = 1;
SS(0) = 1;
for(int i=1; i<N0m; i++){
FF(i) = FF(i-1) + i;
SS(i) = FF(i) + i;}
FF -= 1;
SS -= 1;

arma::rowvec theta(4);
double theta2p;
double theta34;
double theta24;
double theta23;
double A1theta23;
arma::rowvec init = arma::zeros(1,2);
init(0) = A2;
double x1x2;
double B1theta23;
double x1x2theta4;
double u1u2;
double iden;
arma::mat Du(N0,2);
arma::rowvec sapm(2);
arma::rowvec uu(2);

arma::vec mmm(N);
int ip;
int jp;
arma::mat rn1(N0m,2);
arma::mat out = arma::zeros(B,n);
int C1;
int C2;
int E1;
int E2;
int F;
int L;
int N0mm;
N0mm = N0m - 1;

for(int b=0; b<B; b++){

		for(int m=0; m<M; m++){
		
		E1 = m*4;
		E2 = E1 + 3;
		theta = t.submat(b,E1,b,E2);
		theta2p = theta(1) + 1;
		theta24 = theta(1)*theta(3);
		theta23 = theta(1)*theta(2);
		theta34 = theta(2)*theta(3);
		A1theta23 = A1 + theta23;

		init(1) = x(m,1);  
		B1 = x(m,0);
		Du = arma::zeros(N0,2);
		x1x2 = A1 + B1;
		B1theta23 = B1 + theta23;
		x1x2theta4 = theta(3)*x1x2;
		u1u2 = x.at(m,1);
		iden = theta.at(0)/(2*x1x2*u1u2 + theta2p*(x1x2theta4 + theta.at(2)*u1u2) + 2*theta34);
		Du.at(0,0) = (A1*(init.at(1) + theta24) - init.at(0)*B1theta23)*iden;
		Du.at(0,1) = (B1*(init.at(0) + theta24) - init.at(1)*A1theta23)*iden;

			ip = 0;
			//rn1.randn();///////// HERE HERE
			F = m*N0m;
			L = F + N0mm;
			rn1.unsafe_col(0) = LUCY1.submat(F,b,L,b);
			rn1.unsafe_col(1) = LUCY2.submat(F,b,L,b);
			for(int i=0; i<N0m; i++){
			ip += 1;
			sapm = init;
			jp = FF.at(i);
			for(int j=0; j<ip; j++){
			sapm.at(0) += Du.at(j,0)*qrb.at(jp);
			sapm.at(1) += Du.at(j,1)*qrb.at(jp);
			jp += 1;}
			uu = rn1.row(i);
			uu *= ssapv.at(i);
			uu += sapm;
			u1u2 = uu.at(0)+uu.at(1);
			iden = theta.at(0)/(2*x1x2*u1u2 + theta2p*(x1x2theta4 + theta.at(2)*u1u2) + 2*theta34);
			Du.at(ip,0) = (A1*(uu.at(1) + theta24) - uu.at(0)*B1theta23)*iden;
			Du.at(ip,1) = (B1*(uu.at(0) + theta24) - uu.at(1)*A1theta23)*iden;

			}

		mmm = QRBEVAL*Du.col(0);
		mmm += init(0);
		
		C1 = N*m;
		C2 = C1 + Nm;
		out.submat(b,C1,b,C2) = mmm.t();
		
		}

}

return as<NumericMatrix>(wrap(out));

}

RcppExport SEXP nsel4cpp(SEXP THETA, SEXP sig, SEXP mu, SEXP y ) {

Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericVector sigr(sig);
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix yr(y);

int B = THETAr.nrow();
int n = yr.ncol();

arma::mat theta(THETAr.begin(), B, 4, false);
arma::vec SIG(sigr.begin(), B, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);

arma::mat MUt = MU.t();
arma::mat Yt = Y.t();
arma::mat thetat = theta.t();

arma::vec F = arma::zeros(B);
for(int i=0; i<B; i++){
F.at(i) -= dot(MUt.unsafe_col(i),MUt.unsafe_col(i));
F.at(i) /= SIG.at(i);}
F -= n*log(SIG);
F *= 0.5;

arma::vec S = arma::zeros(B);
for(int i=0; i<B; i++){
S.at(i) += dot(Yt.unsafe_col(i),Yt.unsafe_col(i));}
S *= 0.5;

arma::vec ISIG = 1/SIG;

double L;
double U;
arma::mat out = arma::zeros(4,B);
arma::vec Q(B);
arma::vec hy(n);
for(int i=0; i<B; i++){
Q = F;
hy = Yt.unsafe_col(i);
for(int j=0; j<B; j++){
L = dot(hy,MUt.unsafe_col(j));
L -= S.at(i);
L *= ISIG.at(j);
Q.at(j) += L;}
U = max(Q);
Q -= U;
Q = exp(Q);
L = 0;
for(int j=0; j<B; j++){
out.unsafe_col(i) += Q.at(j)*thetat.unsafe_col(j);
L += Q.at(j);}
out.unsafe_col(i) /= L;}

return as<NumericMatrix>(wrap(out.t()));

}


//RcppExport SEXP sigcpp(SEXP sig2, SEXP mu1, SEXP mu2, SEXP y ) {

//Rcpp::NumericVector sig2r(sig2);
//Rcpp::NumericMatrix mu1r(mu1);
//Rcpp::NumericMatrix mu2r(mu2);
//Rcpp::NumericMatrix yr(y);

//int B = yr.nrow();
//int n = yr.ncol();

//arma::vec SIG2(sig2r.begin(), B, false);
//arma::mat MU1(mu1r.begin(), B, n, false);
//arma::mat MU2(mu2r.begin(), B, n, false);
//arma::mat Y(yr.begin(), B, n, false);

//arma::mat MU1t = MU1.t();
//arma::mat MU2t = MU2.t();
//arma::mat Yt = Y.t();

//arma::vec S = -0.5*n*log(SIG2);
//arma::vec ISIG2 = 1/SIG2;

//arma::vec A = arma::zeros(B);
//for(int i=0; i<B; i++){
//A.at(i) = dot(Yt.unsafe_col(i),Yt.unsafe_col(i));}
//A *= 0.5;

//arma::vec BB = arma::zeros(B);
//for(int i=0; i<B; i++){
//BB.at(i) = dot(MU2t.unsafe_col(i),MU2t.unsafe_col(i));}
//BB %= ISIG2;
//BB *= 0.5;

//arma::mat diff = Yt;
//diff -= MU1t;
//diff %= diff;
//arma::vec T = arma::zeros(B);
//for(int i=0; i<B; i++){
//for(int k=0; k<n; k++){
//T.at(i) += diff.at(k,i);}}
//T *= 0.5;

//arma::mat MU2trat = MU2t;
//MU2trat.each_row() %= ISIG2.t();
 
//arma::vec out = arma::zeros(B);
//arma::vec Q1(B);
//arma::vec Q2(B);
//arma::vec hy(n);
//double mK1;
//double mK2;
//double L1;
//double L2;

//for(int i=0; i<B; i++){
//hy = Yt.unsafe_col(i);
//Q1 = S;
//Q2 = S;
//for(int j=0; j<B; j++){
//Q1.at(j) -= ISIG2.at(j)*T.at(i);
//Q2.at(j) -= ISIG2.at(j)*A.at(i);
//Q2.at(j) -= BB.at(j);
//Q2.at(j) += dot(hy,MU2trat.unsafe_col(j));}

//mK1 = max(Q1);
//Q1 -= mK1;
//Q1 = exp(Q1);
//L1 = sum(Q1);
//L1 = log(L1);
//L1 += mK1;
//out.at(i) = L1;
//mK2 = max(Q2);
//Q2 -= mK2;
//Q2 = exp(Q2);
//L2 = sum(Q2);
//L2 = log(L2);
//L2 += mK2;
//out.at(i) -= L2;}

//return as<NumericVector>(wrap(out));

//}






RcppExport SEXP naelcpp(SEXP THETA, SEXP sig, SEXP mu, SEXP y, SEXP order ) {

Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericVector sigr(sig);
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix orderr(order);

int B = THETAr.nrow();
int n = yr.ncol();

arma::mat theta(THETAr.begin(), B, 4, false);
arma::vec SIG(sigr.begin(), B, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);
arma::mat ORD(orderr.begin(), B, 4, false);     

arma::mat MUt = MU.t();
arma::mat Yt = Y.t();

arma::vec F = arma::zeros(B);
for(int i=0; i<B; i++){
F.at(i) -= dot(MUt.unsafe_col(i),MUt.unsafe_col(i));
F.at(i) /= SIG.at(i);}
F -= n*log(SIG);
F *= 0.5;

arma::vec S = arma::zeros(B);
for(int i=0; i<B; i++){
S.at(i) += dot(Yt.unsafe_col(i),Yt.unsafe_col(i));}
S *= 0.5;

arma::umat OT(B,4);
arma::mat theta3 = arma::zeros(B,4);
for(int i=0; i<B; i++){ 
for(int j=0; j<4; j++){ 
OT.at(i,j) = ORD.at(i,j)-1;
theta3.at(i,j) = theta.at(OT.at(i,j),j);}}

double L;
double U;
arma::mat out = arma::zeros(4,B);
arma::vec Q(B);
arma::vec hy(n);
arma::vec L4(B);
double L3;
arma::uvec q1;
arma::uvec q1n(1);
arma::uvec Bm(1);
Bm(0) = B - 1;

arma::vec ISIG = 1/SIG;

for(int i=0; i<B; i++){
Q = F;
hy = Yt.unsafe_col(i);
for(int j=0; j<B; j++){
L = dot(hy,MUt.unsafe_col(j));
L -= S.at(i);
L *= ISIG.at(j);
Q.at(j) += L;}
U = max(Q);
Q -= U;
Q = exp(Q);
L = sum(Q);
Q /= L;

for(int k=0; k<4; k++){
L4 = Q.elem(OT.unsafe_col(k));
L4 = cumsum(L4);
q1 = find(L4<0.5,1,"last");
q1n(0) = q1.n_elem;
	if(q1n(0)>0){
		if(q1(0)<Bm(0)){
		L3 = theta3.at(q1(0),k);
		L3 += theta3.at(q1(0)+1,k);
		L3 *= 0.5;
		} else{
		L3 = theta3.at(Bm(0),k);}
	}
	if(q1n(0)==0){
	L3 = theta3.at(0,k);}
out.at(k,i) = L3;}

}

return as<NumericMatrix>(wrap(out.t()));

}

RcppExport SEXP modcpp(SEXP sig, SEXP mu1, SEXP mu2, SEXP y) {
             
Rcpp::NumericVector sigr(sig);
Rcpp::NumericMatrix mu1r(mu1);   
Rcpp::NumericMatrix mu2r(mu2);
Rcpp::NumericMatrix yr(y);

int B = yr.nrow();
int n = yr.ncol();

arma::vec SIG2(sigr.begin(), B, false);
arma::mat MU1(mu1r.begin(), B, n, false);
arma::mat MU2(mu2r.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);

arma::vec F = -0.5*n*log(SIG2);

arma::mat Yt = Y.t();
arma::mat MU1t = MU1.t();
arma::mat MU2t = MU2.t();

arma::vec MM1(B);
arma::vec MM2(B);
arma::vec YY(B);
for(int i=0; i<B; i++){
YY(i) = dot(Yt.col(i),Yt.col(i));
MM1(i) = dot(MU1t.col(i),MU1t.col(i));
MM2(i) = dot(MU2t.col(i),MU2t.col(i));}
YY *= 0.5;
MM1 *= 0.5;
MM2 *= 0.5;

arma::vec ISIG2 = 1/SIG2;

arma::vec hey = arma::zeros(n);
arma::vec Q1 = arma::zeros(B);
arma::vec Q2 = arma::zeros(B);
double K;
double G1;
double G2;

arma::vec out(B);
out.fill(1.0);

for(int i=0; i<B; i++){

hey = Yt.col(i);
Q1 = F;
Q2 = F;
K = YY(i);

for(int j=0; j<B; j++){

G1 = K;
G1 += MM1(j);
G1 -= dot(hey,MU1t.col(j));
G1 *= ISIG2(j);
Q1(j) -= G1;

G2 = K;
G2 += MM2(j);
G2 -= dot(hey,MU2t.col(j));
G2 *= ISIG2(j);
Q2(j) -= G2;

}

Q1 = exp(Q1);
Q2 = exp(Q2);
G1 = sum(Q1);
G2 = sum(Q2);

if(G1>G2){
out(i) = 1;} else{
out(i) = 2;}

}

return as<NumericVector>(wrap(out));

}

RcppExport SEXP jakstatsolver2cpp(SEXP S, SEXP THETA2, SEXP U0, SEXP TAU, SEXP qrbeval, SEXP QRB, SEXP SAPV, SEXP epo, SEXP randoms1, SEXP randoms2, SEXP randoms3, SEXP randoms4) {

Rcpp::NumericVector Sr(S);
Rcpp::NumericMatrix THETA2r(THETA2);
Rcpp::NumericVector U0r(U0);                 
Rcpp::NumericVector TAUr(TAU);                 
Rcpp::NumericMatrix qrbevalr(qrbeval);
Rcpp::NumericVector QRBr(QRB);
Rcpp::NumericVector SAPVr(SAPV);
Rcpp::NumericVector epor(epo);
Rcpp::NumericMatrix randoms1r(randoms1);
Rcpp::NumericMatrix randoms2r(randoms2);
Rcpp::NumericMatrix randoms3r(randoms3);
Rcpp::NumericMatrix randoms4r(randoms4);

int N0m = SAPVr.size();
int N0 = N0m + 1;
int B = THETA2r.nrow();
int n = qrbevalr.nrow();
int nm = n - 1;
int n2m = 2*n - 1;
int n2 = n*2;
int n3m = n*3 - 1;

arma::vec s(Sr.begin(), N0, false);
arma::mat t(THETA2r.begin(), B, 4, false);
arma::vec u0(U0r.begin(), B, false);
arma::vec tau(TAUr.begin(), B, false);
arma::mat QRBEVAL(qrbevalr.begin(), n, N0, false);
arma::vec qrb(QRBr.begin(), QRBr.size(), false);
arma::vec sapv(SAPVr.begin(), N0m, false);
arma::mat EPO(epor.begin(),B,N0,false); 
arma::mat RA1(randoms1r.begin(),N0m,B,false); 
arma::mat RA2(randoms2r.begin(),N0m,B,false); 
arma::mat RA3(randoms3r.begin(),N0m,B,false); 
arma::mat RA4(randoms4r.begin(),N0m,B,false); 

arma::vec ssapv(N0m);
ssapv = sqrt(sapv);

arma::uvec FF(N0m);
arma::uvec SS(N0m);
FF(0) = 1;
SS(0) = 1;
for(int i=1; i<N0m; i++){
FF(i) = FF(i-1) + i;
SS(i) = FF(i) + i;}
FF -= 1;
SS -= 1;

arma::mat rn1(N0m,4);
arma::rowvec init = arma::zeros(1,4);
arma::mat Du(N0,4);
arma::rowvec sapm1(4);
double delay;
arma::rowvec uu(4);
arma::vec UU = arma::zeros(N0);
arma::vec mmm1(n);
arma::vec mmm2(n);
arma::vec mmm3(n);
int ip;
int jp;
arma::vec sdr(N0);
double msdr;
arma::uvec whee(1);
double tu1;
double tu2;
double tu3;
double tu4;

arma::mat out = arma::zeros(B,3*n);

for(int b=0; b<B; b++){

init(0) = u0(b);
Du = arma::zeros(N0,4);

Du.at(0,0) = -t(b,0)*init(0)*EPO(b,0);
Du.at(0,1) = t(b,0)*init(0)*EPO(b,0);
Du.at(0,2) = 0;
Du.at(0,3) = 0;

			ip = 0;
			//rn1.randn();   												///////// HERE HERE
			rn1.unsafe_col(0) = RA1.unsafe_col(b);
			rn1.unsafe_col(1) = RA2.unsafe_col(b);
			rn1.unsafe_col(2) = RA3.unsafe_col(b);
			rn1.unsafe_col(3) = RA4.unsafe_col(b);	
			for(int i=0; i<N0m; i++){
			ip += 1;
			sapm1 = init;
			jp = FF.at(i);
			for(int j=0; j<ip; j++){
			sapm1.at(0) += Du.at(j,0)*qrb.at(jp);
			sapm1.at(1) += Du.at(j,1)*qrb.at(jp);
			sapm1.at(2) += Du.at(j,2)*qrb.at(jp);
			sapm1.at(3) += Du.at(j,3)*qrb.at(jp);
			jp += 1;}
			uu = rn1.row(i);
			//uu.zeros();
			uu *= ssapv.at(i);
			uu += sapm1;
			tu2 = uu(1);
			tu2 *= tu2;
			tu2 *= t(b,1);
			tu1 = t(b,2)*uu(2);
			tu3 = t(b,0)*uu(0)*EPO(b,ip);
			UU(ip) = uu(3);
			delay = s(ip) - tau(b);
			if(delay<0){
			Du.at(ip,0) = -tu3;
			Du.at(ip,1) = tu3- tu2;
			Du.at(ip,2) = -tu1 + 0.5*tu2;
			Du.at(ip,3) = tu1;
			} else{
			sdr = s - delay;
			sdr %= sdr;
			msdr = min(sdr);
			whee = find(sdr==msdr);
			tu4 = t(b,3)*UU(whee(0));
			Du.at(ip,0) = -tu3 + 2*tu4;
			Du.at(ip,1) = tu3 - tu2;
			Du.at(ip,2) = -tu1 + 0.5*tu2;
			Du.at(ip,3) = tu1 - tu4;
			}
			}
			
mmm1 = QRBEVAL*Du.col(0);
mmm1 += init(0);
mmm2 = QRBEVAL*Du.col(1);
mmm2 += init(1);
mmm3 = QRBEVAL*Du.col(2);
mmm3 += init(2);

out.submat(b,0,b,nm) = mmm1.t();
out.submat(b,n,b,n2m) = mmm2.t();
out.submat(b,n2,b,n3m) = mmm3.t();		
}




return as<NumericMatrix>(wrap(out));

}


RcppExport SEXP jakstatnsel4cpp(SEXP THETA, SEXP mu, SEXP v, SEXP y ) {

Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix yr(y);

int B = THETAr.nrow();
int n = yr.ncol();

arma::mat t(THETAr.begin(), B, 6, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);

arma::vec L = arma::zeros(B);
arma::vec Q = arma::zeros(B);
double U1;
double U2;
for(int b=0; b<B; b++){
for(int k=0; k<n; k++){
L(b) -= log(V(b,k));
U1 = MU(b,k);
U1 *= MU(b,k);
U1 /= V(b,k);
Q(b) -= U1;}}
L *= 0.5;
Q *= 0.5;
arma:: vec LQ = L+Q;

arma::mat Yt = Y.t();
arma::mat F = Yt;
F %= Yt;
F *= 0.5;
arma::mat Bm = 1/V.t();
arma::mat T = MU.t();
T %= Bm;

arma::mat tt = t.t();

arma::vec hy(n);
arma::mat out = arma::zeros(6,B);

for(int i=0; i<B; i++){
hy = Yt.col(i);
U2 = 0;
for(int j=0; j<B; j++){
U1 = LQ(j);
U1 -= dot(F.unsafe_col(i),Bm.unsafe_col(j));
U1 += dot(Yt.unsafe_col(i),T.unsafe_col(j));
U1 = exp(U1);
out.col(i) += U1*tt.col(j);
U2 += U1;}
out.col(i) /= U2;}

return as<NumericMatrix>(wrap(out.t()));

}

RcppExport SEXP jakstatnaelcpp(SEXP THETA, SEXP mu, SEXP v, SEXP y, SEXP order ) {

Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix orderr(order);

int B = THETAr.nrow();
int n = yr.ncol();

arma::mat t(THETAr.begin(), B, 6, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);
arma::mat ORD(orderr.begin(), B, 6, false);     

arma::umat OT(B,6);
arma::mat theta3 = arma::zeros(B,6);
for(int i=0; i<B; i++){ 
for(int j=0; j<6; j++){ 
OT(i,j) = ORD(i,j)-1;
theta3(i,j) = t(OT(i,j),j);}}

arma::vec L = arma::zeros(B);
arma::vec Q = arma::zeros(B);
double U1;
arma::vec U2(B);
arma::vec L4(B);
double L3;
arma::uvec q1;
arma::uvec q1n(1);
arma::uvec Bmr(1);
Bmr(0) = B - 1;
for(int b=0; b<B; b++){
for(int k=0; k<n; k++){
L(b) -= log(V(b,k));
U1 = MU(b,k);
U1 *= MU(b,k);
U1 /= V(b,k);
Q(b) -= U1;}}
L *= 0.5;
Q *= 0.5;
arma:: vec LQ = L+Q;

arma::mat Yt = Y.t();
arma::mat F = Yt;
F %= Yt;
F *= 0.5;
arma::mat Bm = 1/V.t();
arma::mat T = MU.t();
T %= Bm;

arma::mat tt = t.t();

arma::vec hy(n);
arma::mat out = arma::zeros(6,B);

for(int i=0; i<B; i++){
hy = Yt.col(i);

U1 = 0;
for(int j=0; j<B; j++){
U2(j) = LQ(j);
U2(j) -= dot(F.unsafe_col(i),Bm.unsafe_col(j));
U2(j) += dot(Yt.unsafe_col(i),T.unsafe_col(j));
U2(j) = exp(U2(j));
U1 += U2(j);}
U2 /= U1;

for(int k=0; k<6; k++){
L4 = U2.elem(OT.col(k));
L4 = cumsum(L4);
q1 = find(L4<0.5,1,"last");
q1n(0) = q1.n_elem;
	if(q1n(0)>0){
		if(q1(0)<Bmr(0)){
		L3 = theta3(q1(0),k);
		L3 += theta3(q1(0)+1,k);
		L3 *= 0.5;
		} else{
		L3 = theta3(Bmr(0),k);}
	}
	if(q1n(0)==0){
	L3 = theta3(0,k);}
out(k,i) = L3;}

}

return as<NumericMatrix>(wrap(out.t()));

}



RcppExport SEXP jakstatsigcpp(SEXP mu1, SEXP mu2, SEXP v, SEXP y ) {

Rcpp::NumericMatrix mu1r(mu1);                 
Rcpp::NumericMatrix mu2r(mu2);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix yr(y);

int B = yr.nrow();
int n = yr.ncol();

arma::mat MU1(mu1r.begin(), B, n, false);
arma::mat MU2(mu2r.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);

arma::mat MU1t = MU1.t();
arma::mat MU2t = MU2.t();
arma::mat Yt = Y.t();
arma::mat Vt = V.t();

arma::vec F = arma::zeros(B);
for(int j=0; j<B; j++){
for(int k=0; k<n; k++){
F(j) -= log(V(j,k));}}
F *= 0.5;

arma::vec S = arma::zeros(B);
for(int j=0; j<B; j++){
for(int k=0; k<n; k++){
S(j) -= MU2(j,k)*MU2(j,k)/V(j,k);}}
S *= 0.5;
S += F;

arma::mat L = Yt - MU1t;
L %= L;
L *= 0.5;
arma::mat R = Yt;
R %= R;
R *= 0.5;
arma::mat P = 1/Vt;
arma::mat Q = MU2t;
Q %= P;

arma::vec L1;
arma::vec L2;
double K1;
double K2;
double W1;
double W2;
arma::vec out = arma::zeros(B);

for(int i=0; i<B; i++){

L1 = F;
L2 = S;
for(int j=0; j<B; j++){
L1(j) -= dot(L.unsafe_col(i),P.unsafe_col(j));
L2(j) -= dot(R.unsafe_col(i),P.unsafe_col(j));
L2(j) += dot(Q.unsafe_col(j),Yt.unsafe_col(i));}
K1 = max(L1);
K2 = max(L2);
L1 -= K1;
L2 -= K2;
L1 = exp(L1);
L2 = exp(L2);
W1 = sum(L1);
W2 = sum(L2);
W1 = log(W1);
W2 = log(W2);
W1 += K1;
W2 += K2;
out(i) = W1 - W2;}

return as<NumericVector>(wrap(out));

}

RcppExport SEXP fitzsolver2cpp(SEXP S, SEXP THETA2, SEXP qrbeval, SEXP QRB, SEXP SAPV, SEXP randoms1, SEXP randoms2) {

Rcpp::NumericVector Sr(S);
Rcpp::NumericMatrix THETA2r(THETA2); 
Rcpp::NumericMatrix qrbevalr(qrbeval);
Rcpp::NumericVector QRBr(QRB);
Rcpp::NumericVector SAPVr(SAPV);
Rcpp::NumericMatrix randoms1r(randoms1);
Rcpp::NumericMatrix randoms2r(randoms2);

int N0m = SAPVr.size();
int N0 = N0m + 1;
int B = THETA2r.nrow();
int n = qrbevalr.nrow();

arma::vec s(Sr.begin(), N0, false);
arma::mat t(THETA2r.begin(), B, 3, false);
arma::mat QRBEVAL(qrbevalr.begin(), n, N0, false);
arma::vec qrb(QRBr.begin(), QRBr.size(), false);
arma::vec sapv(SAPVr.begin(), N0m, false);
arma::mat Lucy1(randoms1r.begin(), N0m, B, false);
arma::mat Lucy2(randoms2r.begin(), N0m, B, false);

arma::vec ssapv(N0m);
ssapv = sqrt(sapv);

arma::uvec FF(N0m);
arma::uvec SS(N0m);
FF(0) = 1;
SS(0) = 1;
for(int i=1; i<N0m; i++){
FF(i) = FF(i-1) + i;
SS(i) = FF(i) + i;}
FF -= 1;
SS -= 1;

arma::mat Du(N0,2);
arma::rowvec init = arma::zeros(1,2);
init(0) -= 1;
init(1) += 1;
arma::mat rn1(N0m,2);
arma::rowvec sapm(2);
arma::vec mmm1(n);
int ip;
int jp;
arma::rowvec uu(2);
arma::mat out = arma::zeros(B,n);
double V;
double R;
double V3;

for(int b=0; b<B; b++){

Du = arma::zeros(N0,4);
V = init(0);
R = init(1);
V3 = pow(V,3)/3;

Du(0,0) = t(b,2)*(V - V3 + R);
Du(0,1) = -(V - t(b,0) + t(b,1)*R)/t(b,2);

			ip = 0;
			//rn1.randn();   												///////// HERE HERE
			rn1.unsafe_col(0) = Lucy1.unsafe_col(b);
			rn1.unsafe_col(1) = Lucy2.unsafe_col(b);
			for(int i=0; i<N0m; i++){
			ip += 1;
			sapm = init;
			jp = FF.at(i);
			for(int j=0; j<ip; j++){
			sapm.at(0) += Du.at(j,0)*qrb.at(jp);
			sapm.at(1) += Du.at(j,1)*qrb.at(jp);
			jp += 1;}
			uu = rn1.row(i);
			uu *= ssapv.at(i);
			uu += sapm;
			V = uu(0);
			R = uu(1);
			V3 = pow(V,3)/3;

			Du(ip,0) = t(b,2)*(V - V3 + R);
			Du(ip,1) = -(V - t(b,0) + t(b,1)*R)/t(b,2);
		
			}

mmm1 = QRBEVAL*Du.col(0);
mmm1 += init(0);

out.row(b) = mmm1.t();
}

return as<NumericMatrix>(wrap(out));

}


RcppExport SEXP fitzhughnsel4cpp(SEXP THETA, SEXP mu, SEXP v, SEXP y ) {

Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix yr(y);

int B = THETAr.nrow();
int n = yr.ncol();

arma::mat t(THETAr.begin(), B, 3, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);

arma::vec L = arma::zeros(B);
arma::vec Q = arma::zeros(B);
double U1;
double U2;
for(int b=0; b<B; b++){
for(int k=0; k<n; k++){
L(b) -= log(V(b,k));
U1 = MU(b,k);
U1 *= MU(b,k);
U1 /= V(b,k);
Q(b) -= U1;}}
L *= 0.5;
Q *= 0.5;
arma:: vec LQ = L+Q;

arma::mat Yt = Y.t();
arma::mat F = Yt;
F %= Yt;
F *= 0.5;
arma::mat Bm = 1/V.t();
arma::mat T = MU.t();
T %= Bm;

arma::mat tt = t.t();

arma::vec hy(n);
arma::mat out = arma::zeros(3,B);

for(int i=0; i<B; i++){
hy = Yt.col(i);
U2 = 0;
for(int j=0; j<B; j++){
U1 = LQ(j);
U1 -= dot(F.unsafe_col(i),Bm.unsafe_col(j));
U1 += dot(Yt.unsafe_col(i),T.unsafe_col(j));
U1 = exp(U1);
out.col(i) += U1*tt.col(j);
U2 += U1;}
out.col(i) /= U2;}

return as<NumericMatrix>(wrap(out.t()));

}

RcppExport SEXP fitzhughnaelcpp(SEXP THETA, SEXP mu, SEXP v, SEXP y, SEXP order ) {

Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix orderr(order);

int B = THETAr.nrow();
int n = yr.ncol();

arma::mat t(THETAr.begin(), B, 3, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);
arma::mat ORD(orderr.begin(), B, 3, false);     

arma::umat OT(B,3);
arma::mat theta3 = arma::zeros(B,3);
for(int i=0; i<B; i++){ 
for(int j=0; j<3; j++){ 
OT(i,j) = ORD(i,j)-1;
theta3(i,j) = t(OT(i,j),j);}}

arma::vec L = arma::zeros(B);
arma::vec Q = arma::zeros(B);
double U1;
arma::vec U2(B);
arma::vec L4(B);
double L3;
arma::uvec q1;
arma::uvec q1n(1);
arma::uvec Bmr(1);
Bmr(0) = B - 1;
for(int b=0; b<B; b++){
for(int k=0; k<n; k++){
L(b) -= log(V(b,k));
U1 = MU(b,k);
U1 *= MU(b,k);
U1 /= V(b,k);
Q(b) -= U1;}}
L *= 0.5;
Q *= 0.5;
arma:: vec LQ = L+Q;

arma::mat Yt = Y.t();
arma::mat F = Yt;
F %= Yt;
F *= 0.5;
arma::mat Bm = 1/V.t();
arma::mat T = MU.t();
T %= Bm;

arma::mat tt = t.t();

arma::vec hy(n);
arma::mat out = arma::zeros(3,B);

for(int i=0; i<B; i++){
hy = Yt.col(i);

U1 = 0;
for(int j=0; j<B; j++){
U2(j) = LQ(j);
U2(j) -= dot(F.unsafe_col(i),Bm.unsafe_col(j));
U2(j) += dot(Yt.unsafe_col(i),T.unsafe_col(j));
U2(j) = exp(U2(j));
U1 += U2(j);}
U2 /= U1;

for(int k=0; k<3; k++){
L4 = U2.elem(OT.col(k));
L4 = cumsum(L4);
q1 = find(L4<0.5,1,"last");
q1n(0) = q1.n_elem;
	if(q1n(0)>0){
		if(q1(0)<Bmr(0)){
		L3 = theta3(q1(0),k);
		L3 += theta3(q1(0)+1,k);
		L3 *= 0.5;
		} else{
		L3 = theta3(Bmr(0),k);}
	}
	if(q1n(0)==0){
	L3 = theta3(0,k);}
out(k,i) = L3;}

}

return as<NumericMatrix>(wrap(out.t()));

}



RcppExport SEXP fitzhughsigcpp(SEXP mu1, SEXP mu2, SEXP v, SEXP y ) {

Rcpp::NumericMatrix mu1r(mu1);                 
Rcpp::NumericMatrix mu2r(mu2);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix yr(y);

int B = yr.nrow();
int n = yr.ncol();

arma::mat MU1(mu1r.begin(), B, n, false);
arma::mat MU2(mu2r.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);

arma::mat MU1t = MU1.t();
arma::mat MU2t = MU2.t();
arma::mat Yt = Y.t();
arma::mat Vt = V.t();

arma::vec F = arma::zeros(B);
for(int j=0; j<B; j++){
for(int k=0; k<n; k++){
F(j) -= log(V(j,k));}}
F *= 0.5;

arma::vec S = arma::zeros(B);
for(int j=0; j<B; j++){
for(int k=0; k<n; k++){
S(j) -= MU2(j,k)*MU2(j,k)/V(j,k);}}
S *= 0.5;
S += F;

arma::mat L = Yt - MU1t;
L %= L;
L *= 0.5;
arma::mat R = Yt;
R %= R;
R *= 0.5;
arma::mat P = 1/Vt;
arma::mat Q = MU2t;
Q %= P;

arma::vec L1;
arma::vec L2;
double K1;
double K2;
double W1;
double W2;
arma::vec out = arma::zeros(B);

for(int i=0; i<B; i++){

L1 = F;
L2 = S;
for(int j=0; j<B; j++){
L1(j) -= dot(L.unsafe_col(i),P.unsafe_col(j));
L2(j) -= dot(R.unsafe_col(i),P.unsafe_col(j));
L2(j) += dot(Q.unsafe_col(j),Yt.unsafe_col(i));}
K1 = max(L1);
K2 = max(L2);
L1 -= K1;
L2 -= K2;
L1 = exp(L1);
L2 = exp(L2);
W1 = sum(L1);
W2 = sum(L2);
W1 = log(W1);
W2 = log(W2);
W1 += K1;
W2 += K2;
out(i) = W1 - W2;}

return as<NumericVector>(wrap(out));

}

RcppExport SEXP compsolver2cpp(SEXP S, SEXP THETA2, SEXP qrbeval, SEXP QRB, SEXP SAPV, SEXP randoms1, SEXP randoms2) {

Rcpp::NumericVector Sr(S);
Rcpp::NumericMatrix THETA2r(THETA2);
Rcpp::NumericMatrix qrbevalr(qrbeval);
Rcpp::NumericVector QRBr(QRB);
Rcpp::NumericVector SAPVr(SAPV);
Rcpp::NumericMatrix randoms1r(randoms1);
Rcpp::NumericMatrix randoms2r(randoms2);

int N0m = SAPVr.size();
int N0 = N0m + 1;
int B = THETA2r.nrow();
int n = qrbevalr.nrow();

arma::vec s(Sr.begin(), N0, false);
arma::mat t(THETA2r.begin(), B, 3, false);
arma::mat QRBEVAL(qrbevalr.begin(), n, N0, false);
arma::vec qrb(QRBr.begin(), QRBr.size(), false);
arma::vec sapv(SAPVr.begin(), N0m, false);
arma::mat Lucy1(randoms1r.begin(), N0m, B, false);
arma::mat Lucy2(randoms2r.begin(), N0m, B, false);

double D = 400;

arma::vec ssapv(N0m);
ssapv = sqrt(sapv);

arma::uvec FF(N0m);
arma::uvec SS(N0m);
FF(0) = 1;
SS(0) = 1;
for(int i=1; i<N0m; i++){
FF(i) = FF(i-1) + i;
SS(i) = FF(i) + i;}
FF -= 1;
SS -= 1;

arma::mat rn1(N0m,2);
arma::rowvec init = arma::zeros(1,2);
init(0) = D;  
arma::mat Du(N0,2);
arma::rowvec sapm(2);
arma::rowvec uu(2);
arma::mat out = arma::zeros(B,n);
arma::vec mmm(n);
int ip;
int jp;
double rat;

for(int b=0; b<B; b++){

Du = arma::zeros(N0,2);

rat = t(b,1)/t(b,2);
Du.at(0,0) = -t(b,0)*init(0);
Du.at(0,1) = rat*init(0) - t(b,1)*init(1);

			ip = 0;
			//rn1.randn();   												///////// HERE HERE
			rn1.unsafe_col(0) = Lucy1.unsafe_col(b);
			rn1.unsafe_col(1) = Lucy2.unsafe_col(b);
			for(int i=0; i<N0m; i++){
			ip += 1;
			sapm = init;
			jp = FF.at(i);
			for(int j=0; j<ip; j++){
			sapm.at(0) += Du.at(j,0)*qrb.at(jp);
			sapm.at(1) += Du.at(j,1)*qrb.at(jp);
			jp += 1;}
			uu = rn1.row(i);
			uu *= ssapv.at(i);
			uu += sapm;

			Du.at(ip,0) = -t(b,0)*uu(0);
			Du.at(ip,1) = rat*uu(0) - t(b,1)*uu(1);

			}

mmm = QRBEVAL*Du.col(1);
mmm += init(1);
		
out.row(b) = mmm.t();
		
}

return as<NumericMatrix>(wrap(out));

}


RcppExport SEXP compartmentalnsel4cpp(SEXP THETA, SEXP mu, SEXP v, SEXP y ) {

Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix yr(y);

int B = THETAr.nrow();
int n = yr.ncol();

arma::mat t(THETAr.begin(), B, 3, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);

arma::vec L = arma::zeros(B);
arma::vec Q = arma::zeros(B);
double U1;
double U2;
for(int b=0; b<B; b++){
for(int k=0; k<n; k++){
L(b) -= log(V(b,k));
U1 = MU(b,k);
U1 *= MU(b,k);
U1 /= V(b,k);
Q(b) -= U1;}}
L *= 0.5;
Q *= 0.5;
arma:: vec LQ = L+Q;

arma::mat Yt = Y.t();
arma::mat F = Yt;
F %= Yt;
F *= 0.5;
arma::mat Bm = 1/V.t();
arma::mat T = MU.t();
T %= Bm;

arma::mat tt = t.t();

arma::vec hy(n);
arma::mat out = arma::zeros(3,B);

for(int i=0; i<B; i++){
hy = Yt.col(i);
U2 = 0;
for(int j=0; j<B; j++){
U1 = LQ(j);
U1 -= dot(F.col(i),Bm.col(j));
U1 += dot(Yt.col(i),T.col(j));
U1 = exp(U1);
out.col(i) += U1*tt.col(j);
U2 += U1;}
out.col(i) /= U2;}

return as<NumericMatrix>(wrap(out.t()));

}


RcppExport SEXP compartmentalsigcpp(SEXP mu, SEXP v, SEXP y ) {

Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix yr(y);

int B = yr.nrow();
int n = yr.ncol();

arma::mat MU(mur.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);

arma::vec L = arma::zeros(B);
arma::vec Q = arma::zeros(B);
double U1;
double U2;
for(int b=0; b<B; b++){
for(int k=0; k<n; k++){
L(b) -= log(V(b,k));
U1 = MU(b,k);
U1 *= MU(b,k);
U1 /= V(b,k);
Q(b) -= U1;}}
L *= 0.5;
Q *= 0.5;
arma:: vec LQ = L+Q;

arma::mat Yt = Y.t();
arma::mat F = Yt;
F %= Yt;
F *= 0.5;
arma::mat Bm = 1/V.t();
arma::mat T = MU.t();
T %= Bm;

arma::vec hy(n);
arma::vec out = arma::zeros(B);

for(int i=0; i<B; i++){
hy = Yt.col(i);
U2 = 0;
for(int j=0; j<B; j++){
U1 = LQ(j);
U1 -= dot(F.col(i),Bm.col(j));
U1 += dot(Yt.col(i),T.col(j));
U1 = exp(U1);
U2 += U1;}
U2 = log(U2);
out(i) = U2;}

out -= log(B);

return as<NumericVector>(wrap(out));

}


RcppExport SEXP compartmentalnaelcpp(SEXP THETA, SEXP mu, SEXP v, SEXP y, SEXP order ) {

Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix vr(v);
Rcpp::NumericMatrix orderr(order);

int B = THETAr.nrow();
int n = yr.ncol();

arma::mat t(THETAr.begin(), B, 3, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat V(vr.begin(), B, n, false);
arma::mat Y(yr.begin(), B, n, false);
arma::mat ORD(orderr.begin(), B, 4, false);     

arma::umat OT(B,3);
arma::mat theta3 = arma::zeros(B,3);
for(int i=0; i<B; i++){ 
for(int j=0; j<3; j++){ 
OT(i,j) = ORD(i,j)-1;
theta3(i,j) = t(OT(i,j),j);}}

arma::vec L = arma::zeros(B);
arma::vec Q = arma::zeros(B);
double U1;
arma::vec U2(B);
arma::vec L4(B);
double L3;
arma::uvec q1;
arma::uvec q1n(1);
arma::uvec Bmr(1);
Bmr(0) = B - 1;
for(int b=0; b<B; b++){
for(int k=0; k<n; k++){
L(b) -= log(V(b,k));
U1 = MU(b,k);
U1 *= MU(b,k);
U1 /= V(b,k);
Q(b) -= U1;}}
L *= 0.5;
Q *= 0.5;
arma:: vec LQ = L+Q;

arma::mat Yt = Y.t();
arma::mat F = Yt;
F %= Yt;
F *= 0.5;
arma::mat Bm = 1/V.t();
arma::mat T = MU.t();
T %= Bm;

arma::mat tt = t.t();

arma::vec hy(n);
arma::mat out = arma::zeros(3,B);

for(int i=0; i<B; i++){
hy = Yt.col(i);

U1 = 0;
for(int j=0; j<B; j++){
U2(j) = LQ(j);
U2(j) -= dot(F.col(i),Bm.col(j));
U2(j) += dot(Yt.col(i),T.col(j));
U2(j) = exp(U2(j));
U1 += U2(j);}
U2 /= U1;

for(int k=0; k<3; k++){
L4 = U2.elem(OT.col(k));
L4 = cumsum(L4);
q1 = find(L4<0.5,1,"last");
q1n(0) = q1.n_elem;
	if(q1n(0)>0){
		if(q1(0)<Bmr(0)){
		L3 = theta3(q1(0),k);
		L3 += theta3(q1(0)+1,k);
		L3 *= 0.5;
		} else{
		L3 = theta3(Bmr(0),k);}
	}
	if(q1n(0)==0){
	L3 = theta3(0,k);}
out(k,i) = L3;}

}

return as<NumericMatrix>(wrap(out.t()));

}


