 ## Testing the Asymptotic Formula of in the Article:
 ##  Some families of trees arising in permutation analysis
 ##
 ##  Authors: Mathilde Bouvel, Marni Mishna, Cyril Nicaud
 ##  File creation: September 2016

### Read in files
read("CountingSequences.mpl"): ## SList, PList 
read("Approximations.mpl.txt"): ## rho and tau

### Defintions
G:= k->sqrt(rho[k]/2/Pi/Lambda2(tau[k], k)); # gamma[k]
Lambda:= (x,k)->x^2/(1-x) + add(SList[j]*(x/(1-x))^j, j=4..k);
Lambda2:=(x,k)-> subs(X=x, diff(Lambda(X,k),X,X));

### THEOREM 5. For fixed k, the number of prime-degree restricted 
### strong interval trees of size n, denoted P^(k)_n grows asymptotically like:

Pasympt:= (k,n)-> evalf(G(k)/(1-tau[k])^2*rho[k]^(-n)*n^(-3/2));

### Test ratio

for K from 3 to 8 do 
seq( PList[K][n+1]/Pasympt(K,n), n=10..199);
od:

## Test percentage error  (Claim: This is 0.0272821489992)
PList[8][11]-Pasympt(8,10))/PList[8][11];
