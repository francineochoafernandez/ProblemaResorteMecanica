/*Serie de resortes conectados entre sí. Un extremo está fijo a una pared;
mientras que el otro está sujeto a una fuerza constante F. */

/*Se puede emplear un método por elemento finito para determinar los desplazamientos de los resortes
En estos métodos, el dominio de la solución se divide en una malla con puntos discretos o nodos */
//Se aplica la EDP en cada nodo, donde las derivadas parciales se reemplazan por diferencias finitas divididas
/*

Ecuaciones:
F= k x
F1= k1 (x1- x2)
F2= k1 (-x1 + x2) + k2 (x2 -x3)
F3=k2 (-x2 + x3) + k3 (x3 - x4)
F4 = k3 (-x3 + x4) +k4 (x4 -x5)
F5 = k4 (-x4 + x5)

k1          -k1
-k1        k1+k2      -k2
           -k2        k2+k3     -k3
                       -k3       k3+k4      -k4
                                 -k4         k4
*/

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
int i,j;
using namespace std;

typedef struct GaussSeidel //x= [ L+D]^(−1)  * (b− U x)
{
  float a[4][4],br[4][4], x[4],F[4], k1,k2,k3,k4;
  int n, maxit,indi;
  float err;

  void RellenaMat(void)
  {
    for(int i=0; i<4 ; i++)
    {
      for(int j=0; j<4 ; j++)
      {
        a[i][j]=0;
      }
    }

    k1=1;
    k2=1;
    k3=1;
    k4=1;

    a[0][0]= k1 +k2;
    a[0][1]= -k2;

    a[1][0]= -k2;
    a[1][1]= k2+k3;
    a[1][2]= -k3;

    a[2][1]= -k3;
    a[2][2]= k3+k4;
    a[2][3]= -k4;

    a[3][2]= -k4;
    a[3][3]= k4;


    F[0]=0;
    F[1]=0;
    F[2]=0;
    F[3]=1;
    //F[4]=1;

    n=4;
    for(j=0;j<n;j++)
      x[j]=0;


    err=0.001;
    maxit=100;



  }

  void ImprimeOp(void)
  {
    //Imprimiendo matrices
    for(i=0;i<n;i++)
    {
      cout << "|";
      for(j=0;j<n;j++)
      {
        printf("\t%f \t",a[i][j]);
      }
      cout << "|";
      if(i==n/2)
      {
        printf("   *   ");
      }
      else
      {
        printf("       ");
      }

      printf("| %f |\t",x[i]);

      if(i==n/2)
      {
        printf("   =   ");
      }
      else
      {
        printf("       ");
      }

      printf("| %f | \t",F[i]);
      printf("\n");
    }
    printf("\n");
  }

  void MetodoGauss(void)
  {
    int k=0;
    float suma;
    while(k<maxit)
    {
      suma=0;
      k++;
      for(i=0;i<n;i++)
      {
        suma=0;
        for(j=0;j<n;j++)
        {
          if(i!=j)
          {
            suma=suma+a[i][j]*x[j];
          }
          x[i]=(F[i]-suma)/a[i][i];
        }
      }
    }
    printf("\nDespues de %i iteraciones.\n\nLos resultados son:\n",k);
    for(i=0;i<n;i++)
      printf("x[%i]: %f",i+1,x[i]);
      printf(" \n\n " );
  }

  void CopiaMatriz(float mat1[4][4],float mat2[4][4])
  {
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
      {
        mat2[i][j]=mat1[i][j];
      }
    }
  }

  void AcomodandoMatriz()
  {
    float renc, renc2,aux;
    int y,r=0;

    //Combirtiendola en diagonal Dominante
    for(i=0;i<n;i++)
    {
      renc=0;
      for(j=0;j<n;j++)
      {
        if (i!=j)
        {
          renc= renc + abs(br[i][j]);
        }
      }

      y=0;
      r=i;
      if(renc>abs(br[i][i]))//El renglon no cumple así que se cambia orden
      {
        for(int l=0;l<n;l++)//Este for solo es para que se repitan las lineas la n cantidad de veces.
        {//Aquí se hacen todas las posibles combinaciones para ver cual digito es el que iria en la diagonal
          renc2=0;
          aux=0;
          for(int k=0;k<n;k++)
          {
            if(k!=y)
              renc2= renc2 + abs(br[i][k]);
          }

          if(abs(br[i][y])>renc2)
          {
            for(i=0;i<n;i++)
            {
  	           aux=br[i][r];
  	           br[i][r]=br[i][y];
  	           br[i][y]=aux;
            }
          }
          y++;

        }

      }

    }

  }

  void ComparaDiago(int op)
  {
    //Comparando diagonal
    float ren;
    int c=0;

    for(i=0;i<n;i++)
    {
      ren=0;
      for(j=0;j<n;j++)
      {
        if (i!=j)
        {
          ren= ren + abs(br[i][j]);
        }
      }

      if(ren>abs(br[i][i]))
      {
        if(op==1)
        {
          printf("\nLa matriz NO es dominante, no va a convergir.\n");
          printf("\nSe tratara de mover la matriz.\n");
          AcomodandoMatriz();
        }
        else
        {
          printf("\nLa matriz NO se pudo acomodar para volverla dominante.\n");
          indi=3;
        }
        i=n;
        j=n;
        c++;
      }

    }
    if(c==0)
    {
      printf("\nLa matriz SI es dominante, sí va a convergir.\n");
      indi=1;
      if(op==2)
      {
        CopiaMatriz( br, a);
        ImprimeOp();
      }

    }
  }

}GS;

int main(void)
{
  GS sisecs;
  sisecs.RellenaMat();
  sisecs.ImprimeOp();
  sisecs.CopiaMatriz( sisecs.a, sisecs.br);
  sisecs.ComparaDiago(1);
  if(sisecs.indi!=1)
  {
    sisecs.ComparaDiago(2);
  }
  sisecs.MetodoGauss();

  return 0;
}
