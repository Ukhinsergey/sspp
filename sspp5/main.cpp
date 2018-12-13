#include "mpi.h"
#include <stdint.h>
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;


int main(int argc, char **argv) {
	if (argc != 4) {
		cout << "A.dat B.dat C.dat" << endl;
		return 0;
	}

	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	int tagmatra = 1;
	int tagmatrb = 2;
	int tagsizebuf = 3;
	
	int newsize = pow(size, 1.0001/3);
	MPI_Comm cube, square, linex, liney, linez;
	int dims[3] = {newsize, newsize, newsize};
 	int period[3] = {0,0,0}; 
 	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, false, &cube);
 	MPI_Comm_rank(cube, &rank);

 	int coords[3];
 	MPI_Cart_coords(cube, rank, 3, coords);
 	int remain_dims[3] = {1,1,0};
 	MPI_Cart_sub(cube, remain_dims, &square);
 	remain_dims[1] = 0;
 	MPI_Cart_sub(cube, remain_dims, &linex);
 	remain_dims[0] = 0;
 	remain_dims[1] = 1;
 	MPI_Cart_sub(cube, remain_dims, &liney);
 	remain_dims[1] = 0;
 	remain_dims[2] = 1;
 	MPI_Cart_sub(cube, remain_dims, &linez);
 	int sizex, sizey, sizez;
 	MPI_Comm_size(linex, &sizex);
 	MPI_Comm_size(liney, &sizey);
 	MPI_Comm_size(linez, &sizez);
 	MPI_Datatype filetype;
 	uint64_t ma, na, mb, nb;

 	if (coords[2] == 0) {
 		MPI_File file;
 		MPI_File_open(square, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
		char type;
 		MPI_File_read_all(file, &type, 1, MPI_CHAR, &status);
 		MPI_File_read_all(file, &ma, 1, MPI_UNSIGNED_LONG_LONG, &status);
 		MPI_File_read_all(file, &na, 1, MPI_UNSIGNED_LONG_LONG, &status);
 		int n = na;
 		int m = ma;
 		int partx;
 		int party;
 		if ( coords[0] != sizex - 1) {
 			partx = na / sizex;
 		} else {
 			partx = na / sizex + na % sizex;
 		}

 		if ( coords[1] != sizey - 1) {
 			party = ma / sizey;
 		} else {
 			party = ma / sizey + ma % sizey;
 		}

 		int startx = coords[0] * na / sizex;
 		MPI_Type_create_subarray(1,  &n, &partx, &startx, MPI_ORDER_C, MPI_DOUBLE, &filetype);
 		MPI_Type_commit(&filetype);

 		int temp = ma /sizey;
 		int offset = 2 * sizeof(uint64_t) + sizeof (char) + sizeof(double) * n * temp  * coords[1] ; 
 		MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
 		double matra[party][partx];
 		for(int i = 0 ; i < party; ++i) {
 			MPI_File_read(file, &matra[i][0], partx, MPI_DOUBLE, &status);
 		}
 		MPI_File_close(&file);
 		
 		int sendtorank;
 		int sendtocoord[3] = {coords[0], coords[1], coords[1]};
 		MPI_Cart_rank(cube, sendtocoord, &sendtorank);
 		int bufsizes[2] = {partx, party};
 		MPI_Send(bufsizes, 2, MPI_INT, sendtorank, tagsizebuf, cube);
 		for(int j = 0 ; j < party; ++j) {
 			MPI_Send(&matra[j][0], partx, MPI_DOUBLE, sendtorank, tagmatra, cube);
 		}
 			 
 	}
 	int bufsizes[2] = { -1 , -1 };
 	int partx, party;

 	if (coords[1] == coords[2]) {
 		int sendfrom;
 		int sendfromcoords[3] = {coords[0],coords[1],0};
 		MPI_Cart_rank(cube, sendfromcoords, &sendfrom);
 		MPI_Recv(&bufsizes, 2, MPI_INT, sendfrom, tagsizebuf, cube, &status);
 		partx = bufsizes[0];
 		party = bufsizes[1];
 	}

 	int bcastroot;
 	int bcastcoords[1] = {coords[2]} ;
 	MPI_Cart_rank(liney, bcastcoords, &bcastroot);
 	MPI_Bcast(bufsizes, 2, MPI_INT, bcastroot, liney);
 	partx = bufsizes[0];
 	party = bufsizes[1];
 	double matra[party][partx];

 	if (coords[1] == coords[2]) {

 		int sendfrom;
 		int sendfromcoords[3] = {coords[0],coords[1],0};
 		MPI_Cart_rank(cube, sendfromcoords, &sendfrom);

 		for(int j = 0 ; j < party; ++j) {
 			MPI_Recv(&matra[j][0], partx, MPI_DOUBLE, sendfrom, tagmatra, cube, &status);
 		}
 	}

 	for(int j = 0 ; j < party; ++j) {
 		MPI_Bcast(&matra[j][0], partx, MPI_DOUBLE, bcastroot, liney);
 	}


 	//matr A ready


 	if ( coords[2] == 0) {
 		MPI_File file;
 		MPI_File_open(square, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
		char type;
 		MPI_File_read_all(file, &type, 1, MPI_CHAR, &status);
 		MPI_File_read_all(file, &mb, 1, MPI_UNSIGNED_LONG_LONG, &status);
 		MPI_File_read_all(file, &nb, 1, MPI_UNSIGNED_LONG_LONG, &status);
 		int n = nb;
 		int m = mb;
 		int bpartx;
 		int bparty;
 		if ( coords[0] != sizex - 1) {
 			bpartx = nb / sizex;
 		} else {
 			bpartx = nb / sizex + nb % sizex;
 		}

 		if ( coords[1] != sizey - 1) {
 			bparty = mb / sizey;
 		} else {
 			bparty = mb / sizey + mb % sizey;
 		}

 		int startx = coords[0] * nb / sizex;
 		MPI_Type_create_subarray(1,  &n, &bpartx, &startx, MPI_ORDER_C, MPI_DOUBLE, &filetype);
 		MPI_Type_commit(&filetype);

 		int temp = mb /sizey;
 		int offset = 2 * sizeof(uint64_t) + sizeof (char) + sizeof(double) * n * temp  * coords[1] ; 
 		MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
 		double matrb[bparty][bpartx];
 		for(int i = 0 ; i < bparty; ++i) {
 			MPI_File_read(file, &matrb[i][0], bpartx, MPI_DOUBLE, &status);
 		}
 		MPI_File_close(&file);

 		/*if (coords[0] == sizex - 2 && coords[1] == sizey - 1) {
 			for(int i = 0 ; i < party; ++i) {
 				for (int j = 0 ; j < partx; ++j) {
 					cout << matrb[i][j] << ' ';
 				}
 				cout << endl;
 			}
 		}*/

 		int sendtorank;
 		int sendtocoord[3] = {coords[0], coords[1], coords[0]};
 		MPI_Cart_rank(cube, sendtocoord, &sendtorank);
 		int bbufsizes[2] = {bpartx, bparty};
 		MPI_Send(bbufsizes, 2, MPI_INT, sendtorank, tagsizebuf, cube);
 		for(int j = 0 ; j < bparty; ++j) {
 			MPI_Send(&matrb[j][0], bpartx, MPI_DOUBLE, sendtorank, tagmatrb, cube);
 		}
 	}

 	int bbufsizes[2] = { -1 , -1 };
 	int bpartx, bparty;

 	if (coords[0] == coords[2]) {

 		int sendfrom;
 		int sendfromcoords[3] = {coords[0],coords[1],0};
 		MPI_Cart_rank(cube, sendfromcoords, &sendfrom);
 		MPI_Recv(&bbufsizes, 2, MPI_INT, sendfrom, tagsizebuf, cube, &status);
 		bpartx = bbufsizes[0];
 		bparty = bbufsizes[1];
 	}

 	MPI_Cart_rank(linex, bcastcoords, &bcastroot);
 	MPI_Bcast(bbufsizes, 2, MPI_INT, bcastroot, linex);
 	bpartx = bbufsizes[0];
 	bparty = bbufsizes[1];
 	double matrb[bparty][bpartx];

 	if (coords[0] == coords[2]) {

 		int sendfrom;
 		int sendfromcoords[3] = {coords[0],coords[1],0};
 		MPI_Cart_rank(cube, sendfromcoords, &sendfrom);

 		for(int j = 0 ; j < bparty; ++j) {
 			MPI_Recv(&matrb[j][0], bpartx, MPI_DOUBLE, sendfrom, tagmatrb, cube, &status);
 		}
 	}

 	for(int j = 0 ; j < bparty; ++j) {
 		MPI_Bcast(&matrb[j][0], bpartx, MPI_DOUBLE, bcastroot, linex);
 	}


 	/*if (coords[0] == sizex - 2 && coords[1] == sizey - 1 && coords[2] == sizez - 2) {
 		for(int i = 0 ; i < bparty; ++i) {
 			for (int j = 0 ; j < bpartx; ++j) {
 				cout << matrb[i][j] << ' ';
 			}
 			cout << endl;
 		}
 	}*/
 	
 	//matr b ready

 	/*double matrclocal[party][bpartx];
 	if(partx != bparty) {
 		//cout << "cant multyply " << partx << ' ' <<bparty << ' ' << coords[0] << ' ' << coords[1] << ' ' <<coords[2] <<endl;
 		//exit(0);
 	} else {
 		
 		for(int i = 0 ; i < party; ++i) {
 			for(int j = 0 ; j < bpartx; ++j) {
 				matrclocal[i][j] = 0;
 			}
 		}
 		for(int i = 0 ; i < party; ++i) {
 			for(int k = 0 ; k < bparty; ++k) {
 				for (int j = 0; j < bpartx; ++j) {
 					matrclocal[i][j] += matra[i][k] * matrb[k][j];
 				}
 			}
 		}
 	}

 	if (coords[0] == sizex - 2 && coords[1] == sizey - 1 && coords[2] == sizez - 2) {
 		for(int i = 0 ; i < party; ++i) {
 			for (int j = 0 ; j < partx; ++j) {
 				cout << matra[i][j] << ' ';
 			}
 			cout << endl;
 		}
 		cout << endl << "B : " <<endl;
 	}
 	if (coords[0] == sizex - 2 && coords[1] == sizey - 1 && coords[2] == sizez - 2) {
 		for(int i = 0 ; i < bparty; ++i) {
 			for (int j = 0 ; j < bpartx; ++j) {
 				cout << matrb[i][j] << ' ';
 			}
 			cout << endl;
 		}
 		cout << endl << "C : " <<endl;
 	}
 	if (coords[0] == sizex - 2 && coords[1] == sizey - 1 && coords[2] == sizez - 2) {
 		for(int i = 0 ; i < party; ++i) {
 			for (int j = 0 ; j < bpartx; ++j) {
 				cout << matrclocal[i][j] << ' ';
 			}
 			cout << endl;
 		}
 		cout << partx << bparty << endl;
 	}
 	*/
	MPI_Comm_free(&cube);
	MPI_Comm_free(&square);
	MPI_Comm_free(&linez);
	MPI_Comm_free(&linex);
	MPI_Comm_free(&liney);
	MPI_Finalize();
	return 0;
}
