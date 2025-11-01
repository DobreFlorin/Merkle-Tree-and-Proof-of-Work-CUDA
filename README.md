Explicatii blas, neopt, opt_m:

neopt: varianta neoptimizata a calcului aloca memorie  pt matricile N*N , realizeaza operatiile si apoi elibereaza memoria:
blas:  cu ajutorul functiei cblas_dgemm se calculeaza C = B × Aᵗ si D = Cᵗ × A. Bucla se optmizeaza folosind cblas_dgemv care inlocuieste buclele for 
opt_m: codul este optimizat in modul uramtor: acces imbunatatit la memorie(se precalculeaza transpusele), folosirea buclelor in mod optim din p.d.v al ordinii indicilor(B[i*N + k] este constant în bucla internă => scoatere în variabilă (B_ik)
At[k*N + j] are acces secvențial în j, îmbunătățind localitatea memoriei), evitarea recomputarilor inutile; folosirea unui buffer(x_work) in locul lui x.


Explicatii despre efectul acestor optimizari:
Efectul optimizarilor conduce la  un timp mai rapid cu aproape 25% in cazul variantei optimizate.

Explicatii pentru valorile obinute cu cachegrind((I refs, D refs, Branches)

Se observa in fisierele .cache ca in cazul variantei optimizate nr de 
instructiuni( I refs), accesul la date( D refs), d1 cache misses este mult 
mai mic(in cazul d1 cache numarul este cu aproape 80 % mai mic ceea ce innseamna ca versiunea optimizata:
folosește mai bine cache-ul, deoarece datele sunt accesate secvențial 
transpunerile ajută ca liniile cache să fie utilizate eficient, iar in cazul I refs si D refs scaderea se
datoreaza buclelor si precomputarii operatiilor) 


Explicatii pentru existenta buclei de dimensiune N:
Probabil petru a transforma x in forma dorita inainte de calculul final


Graficul de asemena arata  impactul optimizarilor:

-neopt: crește rapid cu cand crste N – versiunea brută, fără optimizări.
-opt_m: mai eficientă, reduce semnificativ timpul, dar încă folosește calcule cu 
-blas: cea mai performanta, de zeci de ori mai rapida, deoarece foloseste BLAS pentru matrici si vectori(graficul este aproape asimptotic)
