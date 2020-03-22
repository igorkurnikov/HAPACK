#include <stringc.h>
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <numer.h>

main()
{
/*
 Memory_Manager Heap(500);
 memory_manager=&Heap;
 Heap.on();
 Heap.mark("testmat");
// Heap.intercept(6);
*/
 {
 int i,j;
 Matrix A;
 Matrix B(5,6,1,4), C(4,1,4,1);
 for(i=1; i<=5; i++)
   for(j=4; j<=9; j++)
     {
       B(i,j)=(j-4)+6*(i-1);
     }
 for(i=4; i<=7; i++)
     {
       C(i,1)=i;
     }
 // Heap.off();
 cout<<"Ausgabe von A, B und C \n"<<A<<B<<C;
 cin>>i;
 // Heap.on();

 cin>>i;
 A.set(4);
 Matrix D(4,3);
 D.set(5);
 cout<<"A.set(4) \n"<<A<<"D.set(5) \n"<<D;
 B.rebase(1,2);
 cout<<"B.rebase(1,2) \n"<<B;
 cout<<"B.get(1,4): "<<B.get(1,4)<<"\n";
 cout<<"B.get(1,2): "<<B.get(1,2)<<"\n";
 cin>>i;
 cout<<"Offset von A off1 off2: "<<A.offset1()<<"  "<<A.offset2()<<"\n";
 cout<<"Offset von B off1 off2: "<<B.offset1()<<"  "<<B.offset2()<<"\n";
 cout<<"Size von A size1 size2: "<<A.size1()<<"  "<<A.size2()<<"\n";
 cout<<"Size von B size1 size2: "<<B.size1()<<"  "<<B.size2()<<"\n";
 cin>>i;
 A=B;
 cout<<"A=B: \n"<<A;
 Matrix E=C;
 cout<<"Matrix E=C: \n"<<E;
 cin>>i;
 A+=A;
 cout<<"A+=A: \n"<<A;
 E-=E;
 cout<<"E-=E: \n"<<E;
 cin>>i;
 B-=2;
 cout<<"B-=2: \n"<<B;
 C+=3;
 cout<<"C+=3: \n"<<C;
 cin>>i;
 B*=3;
 cout<<"B*=3: \n"<<B;
 C/=2.0;
 cout<<"C/=2: \n"<<C;
 cin>>i;
 B=A+B;
 cout<<"B=A+B: \n"<<B;
 C=C+C;
 cout<<"C=C+C: \n"<<C;
 cin>>i;
 // Heap.mark("transpose");
 A.reset(4,4,1,2);
 for(i=1; i<5; i++)
   for(j=2; j<6; j++)
     A(i,j)=(i-1)*4+(j-2);
 cout<<"A und A.transpose(): \n"<<A;
 A.transpose();
 cout<<A;
// Heap.print(cout);
// Heap.print_memory_usage(cout);
 // Heap.mark("testmat1");
 cin>>i;
 B.reset(4,3,0,3);
 for(i=0; i<4; i++)
   for(j=3; j<6; j++)
     B(i,j)=(i)*3+(j-3);
 cout<<B;
 C.reset(4,3,4,1);
 C.set(0);
 C.multiply(A,B);
 cout<<"C.muliply(A,B):\n"<<C;
 Matrix F=A*B;
 cout<<"F=A*B: \n"<<F;
 A.write("matout");
 B.write("matout");
 A.set(0);
 B.set(0);
// A.read("matout");
 B.read("matout");
 cout<<"Schreiben und Lesen in bzw. aus einem File,\n \
	(A und B zwischenzeitlich auf Null gesetzt): \n"<<A<<B;
 cin>>i;
 fstream file("matout",ios::in | ios::out);
 B.write(file);
 file.close();
 B.set(0);
 file.open("matout",ios::in | ios::out);
 B.read(file);
 file.close();
 cout<<"Schreiben und Lesen in bzw. aus einem File,\n \
	B (zwischenzeitlich auf Null gesetzt): \n"<<B;
 cin>>i;
 cout<<"C und C.transpose(): \n"<<C;
 C.transpose();
 cout<<C;
 cout<<"B und B.transpose(): \n"<<B;
 B.transpose();
 cout<<B;
 cin>>i;
 }
 // Heap.print(cout);
 return 0;
}
