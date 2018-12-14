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


	double timefile = 0;
	double timetemp1, timetemp2;
	double realtime = 0;
	double realtimetemp1, realtimetemp2;

    double *matra, *matrb;
    int partx, party;
    int bpartx, bparty;

	if (coords[2] == 0) {
		timetemp1 = MPI_Wtime();
		MPI_File file;
		MPI_File_open(square, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
		char type;
		MPI_File_read_all(file, &type, 1, MPI_CHAR, &status);
		MPI_File_read_all(file, &ma, 1, MPI_UNSIGNED_LONG_LONG, &status);
		MPI_File_read_all(file, &na, 1, MPI_UNSIGNED_LONG_LONG, &status);
		timetemp2 = MPI_Wtime();
		timefile += timetemp2 - timetemp1;

		realtimetemp1 = MPI_Wtime();

		int n = na;
		int m = ma;
		if ( coords[0] != sizex - 1) {
			partx = ma / sizex;
		} else {
			partx = ma / sizex + ma % sizex;
		}

		if ( coords[1] != sizey - 1) {
			party = na / sizey;
		} else {
			party = na / sizey + na % sizey;
		}

		int starty = coords[1] * (na / sizey);
		MPI_Type_create_subarray(1,  &n, &party, &starty, MPI_ORDER_C, MPI_DOUBLE, &filetype);
		MPI_Type_commit(&filetype);

		int offset = 2 * sizeof(uint64_t) + sizeof (char) + sizeof(double) * coords[0] * na * (ma / sizex); 

		realtimetemp2 = MPI_Wtime();
		realtime += realtimetemp2 - realtimetemp1;

		timetemp1 = MPI_Wtime();
		MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);


		matra = new double[partx * party];


		MPI_File_read(file, matra , partx * party, MPI_DOUBLE, &status);

		MPI_File_close(&file);

		timetemp2 = MPI_Wtime();
		timefile += timetemp2 - timetemp1;


		realtimetemp1 = MPI_Wtime();

		int sendtorank;
		int sendtocoord[3] = {coords[0], coords[1], coords[1]};
		MPI_Cart_rank(cube, sendtocoord, &sendtorank);
		int bufsizes[2] = {partx, party};
		MPI_Send(bufsizes, 2, MPI_INT, sendtorank, tagsizebuf, cube);
        if( coords[1] != 0) {
            MPI_Send(matra, partx * party, MPI_DOUBLE, sendtorank, tagmatra, cube); 
        }


        realtimetemp2 = MPI_Wtime();
        realtime = realtimetemp2 - realtimetemp1;


    }

    realtimetemp1 = MPI_Wtime();

    int bufsizes[2] = { -1 , -1 };

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

    if(coords[1] == 0 && coords[2] == 0) {
    } else {
        matra = new double[partx * party];
    }


    if (coords[1] == coords[2]) {

        int sendfrom;
        int sendfromcoords[3] = {coords[0],coords[1],0};
        MPI_Cart_rank(cube, sendfromcoords, &sendfrom);
        if(coords[2] != 0) {
            MPI_Recv(matra, partx * party, MPI_DOUBLE, sendfrom, tagmatra, cube, &status);
        }

    }

    MPI_Bcast(matra, partx * party, MPI_DOUBLE, bcastroot, liney);

    /*if(coords[0] == 2 && coords[1] == 2 && coords[2] == 2) {
        for(int i = 0 ; i < partx; ++i) {
            for (int j = 0 ; j < party ; ++j) {
                cout << matra[i * party + j] << ' ';
            }
            cout << endl;
        }

    }*/
 	//matr A ready

    realtimetemp2 = MPI_Wtime();
    realtime = realtimetemp2 - realtimetemp1;

    if ( coords[2] == 0) {
        timetemp1 = MPI_Wtime();
        MPI_File file;
        MPI_File_open(square, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
        char type;
        MPI_File_read_all(file, &type, 1, MPI_CHAR, &status);
        MPI_File_read_all(file, &mb, 1, MPI_UNSIGNED_LONG_LONG, &status);
        MPI_File_read_all(file, &nb, 1, MPI_UNSIGNED_LONG_LONG, &status);
        timetemp2 = MPI_Wtime();
        timefile += timetemp2 - timetemp1;

        realtimetemp1 = MPI_Wtime();

        int n = nb;
        int m = mb;
        if ( coords[0] != sizex - 1) {
            bpartx = mb / sizex;
        } else {
            bpartx = mb / sizex + mb % sizex;
        }

        if ( coords[1] != sizey - 1) {
            bparty = nb / sizey;
        } else {
            bparty = nb / sizey + nb % sizey;
        }

        int starty = coords[1] * nb / sizey;
        MPI_Type_create_subarray(1,  &n, &bparty, &starty, MPI_ORDER_C, MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);

        int offset = 2 * sizeof(uint64_t) + sizeof (char) + sizeof(double) * n *(mb / sizex)  * coords[0] ; 

        realtimetemp2 = MPI_Wtime();
        realtime = realtimetemp2 - realtimetemp1;

        timetemp1 = MPI_Wtime();
        MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

        matrb = new double[bpartx * bparty];

        MPI_File_read(file, matrb, bpartx * bparty, MPI_DOUBLE, &status);

        MPI_File_close(&file);

        /*if(coords[0] == 2 && coords[1] == 2) {
            for(int i = 0 ; i < bpartx; ++i) {
                for (int j = 0 ; j < bparty ; ++j) {
                    cout << matrb[i * bparty + j] << ' ';
                }
                cout << endl;
            }
        }*/
        timetemp2 = MPI_Wtime();
        timefile = timetemp2 - timetemp1;


        int sendtorank;
        int sendtocoord[3] = {coords[0], coords[1], coords[0]};
        MPI_Cart_rank(cube, sendtocoord, &sendtorank);
        int bbufsizes[2] = {bpartx, bparty};
        MPI_Send(bbufsizes, 2, MPI_INT, sendtorank, tagsizebuf, cube);
        if( coords[0] != 0) {
            MPI_Send(matrb, bpartx * bparty, MPI_DOUBLE, sendtorank, tagmatrb, cube); 
        }

    }

    realtimetemp1 = MPI_Wtime();

    int bbufsizes[2] = { -1 , -1 };

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
    if(coords[0] == 0 && coords[2] == 0) {
    } else {
        matrb = new double [bpartx * bparty];     
    }
    if (coords[0] == coords[2]) {

        int sendfrom;
        int sendfromcoords[3] = {coords[0],coords[1],0};
        MPI_Cart_rank(cube, sendfromcoords, &sendfrom);
        if(coords[2] != 0) {
            MPI_Recv(matrb, bpartx * bparty, MPI_DOUBLE, sendfrom, tagmatrb, cube, &status);
        }
    }

    MPI_Bcast(matrb, bpartx * bparty, MPI_DOUBLE, bcastroot, linex);

    /*if(coords[0] == 0 && coords[1] == 2 && coords[2] == 2) {
        for(int i = 0 ; i < bpartx; ++i) {
            for (int j = 0 ; j < bparty ; ++j) {
                cout << matrb[i * bparty + j] << ' ';
            }
            cout << endl;
        }
    }*/

 	// matr b reade

 
    double *matrclocal = new double [partx * bparty];



    if(bpartx != party) {
        cout << "cant multyply " << bpartx << ' ' <<party << ' ' << coords[0] << ' ' << coords[1] << ' ' <<coords[2] <<endl;
     		//exit(0);
    } else {
        for(int i = 0 ; i < partx; ++i) {
            for(int j = 0 ; j < bparty; ++j) {
                matrclocal[i * bparty + j] = 0;
            }
        }
        for(int i = 0 ; i < partx; ++i) {
            for(int k = 0 ; k < bpartx; ++k) {
                for (int j = 0; j < bparty; ++j) {
                    matrclocal[i * bparty + j] += matra[i * party + k] * matrb[k * bparty + j];
                }
            }
        }
    }




    double *matrc = new double [partx * bparty];

    int rankcmatr;
    int rootcoord[1] = { 0 } ;
    MPI_Cart_rank(linez, rootcoord, &rankcmatr);


    MPI_Reduce(matrclocal, matrc, partx * bparty, MPI_DOUBLE, MPI_SUM, rankcmatr, linez);


    MPI_Comm_rank(cube, &rank);

    realtimetemp2 = MPI_Wtime();
    realtime = realtimetemp2 - realtimetemp1;

    if (rank == 0) {
        timetemp1 = MPI_Wtime();

        MPI_File file;
        MPI_File_open(MPI_COMM_SELF, argv[3], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
        char type = 'd';
        MPI_File_write(file, &type, 1, MPI_CHAR, &status);
        MPI_File_write(file, &ma, 1, MPI_UNSIGNED_LONG_LONG, &status);
        MPI_File_write(file, &nb, 1, MPI_UNSIGNED_LONG_LONG, &status);
        MPI_File_close(&file);	

        timetemp2 = MPI_Wtime();
        timefile = timetemp2 - timetemp1;
    }
    if (coords[2] == 0) {

        MPI_File file;
        MPI_File_open(square, argv[3], MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
        if ( coords[0] != sizex - 1) {
            partx = ma / sizex;
        } else {
            partx = ma / sizex + ma % sizex;
        }

        if ( coords[1] != sizey - 1) {
            party = nb / sizey;	
        } else {
            party = nb / sizey + nb % sizey;
        }
        int n = nb;
        int starty = coords[1] * (nb / sizey);
        MPI_Type_create_subarray(1,  &n, &party, &starty, MPI_ORDER_C, MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);

        int offset = 2 * sizeof(uint64_t) + sizeof (char) + sizeof(double) * n * (ma / sizex) * coords[0]; 


        timetemp1 = MPI_Wtime();
        MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

        MPI_File_write(file, matrc, partx * bparty, MPI_DOUBLE, &status);

        MPI_File_close(&file);
        timetemp2 = MPI_Wtime();
        timefile = timetemp2 - timetemp1;
    }

    double maxtimefile;
    double maxrealtime;
    MPI_Reduce(&timefile, &maxtimefile, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&realtime, &maxrealtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if( rank == 0 ) {
        cout << "maxtimefile = " << maxtimefile << endl << "maxrealtime = " << maxrealtime << endl;
    }


    delete []matrclocal;
    delete []matra;
    delete []matrb;
    MPI_Comm_free(&cube);
    MPI_Comm_free(&square);
    MPI_Comm_free(&linez);
    MPI_Comm_free(&linex);
    MPI_Comm_free(&liney);
    MPI_Finalize();
    return 0;
}
