/*

/*
 *
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 *
 */


#include "stencil_template_parallel.h"

// --- GLOBAL VARIABLE DEFINITIONS ---
// This is where the global variables are actually created.
int g_n_omp_threads = 1;
double* g_per_thread_comp_time = NULL;
__thread double thread_local_comp_time = 0.0;

// ------------------------------------------------------------------
// ------------------------------------------------------------------

int main(int argc, char **argv)
{
  MPI_Comm myCOMM_WORLD;
  int  Rank, Ntasks;
  uint neighbours[4];

  int  Niterations;
  int  periodic;
  vec2_t S, N;
  
  int      Nsources;
  int      Nsources_local;
  vec2_t  *Sources_local;
  double   energy_per_source;

  plane_t   planes[2];  
  buffers_t buffers[2];
  
  int output_energy_stat_perstep;
  
  /* initialize MPI envrionment */
  {
    int level_obtained;
    
    // NOTE: change MPI_FUNNELED if appropriate
    //
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
    if ( level_obtained < MPI_THREAD_FUNNELED ) {
      printf("MPI_thread level obtained is %d instead of %d\n",
	     level_obtained, MPI_THREAD_FUNNELED );
      MPI_Finalize();
      exit(1); }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
    MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD);
  }
  
    // --- SETUP FOR PER-THREAD TIMING ---
    // Get the number of threads that will be used.
    #pragma omp parallel
    {
        #pragma omp master
        { g_n_omp_threads = omp_get_num_threads(); }
    }
    // Allocate memory for the global timing array and initialize to zero.
    g_per_thread_comp_time = (double*)calloc(g_n_omp_threads, sizeof(double));


  /* argument checking and setting */
  int ret = initialize ( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,
			 neighbours, &Niterations,
			 &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
			 &planes[0], &buffers[0] );


  if ( ret )
    {
      printf("task %d is opting out with termination code %d\n",
	     Rank, ret );
      
      MPI_Finalize();
      return 0;
    }
  
  
  int current = OLD;
  double t_start, t_end;
  double total_comm_time = 0.0;
  double total_comp_time = 0.0;

  
  t_start = MPI_Wtime();


for (int iter = 0; iter < Niterations; ++iter)
    
    {
      double section_start_time;
 
      
      /* new energy from sources */
      inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N );
      
      section_start_time = MPI_Wtime();

      /* -------------------------------------- */

      double* current_plane = planes[current].data;
      const int sizex = planes[current].size[_x_];
      const int sizey = planes[current].size[_y_];
      const int full_sizex = sizex + 2;

      MPI_Request requests[8];
      int request_count = 0;

      // [A] Pack East/West buffers
      if (neighbours[WEST] != MPI_PROC_NULL) {
          for (int i = 0; i < sizey; i++) buffers[SEND][WEST][i] = current_plane[(i + 1) * full_sizex + 1];
      }
      if (neighbours[EAST] != MPI_PROC_NULL) {
          for (int i = 0; i < sizey; i++) buffers[SEND][EAST][i] = current_plane[(i + 1) * full_sizex + sizex];
      }


      // [B] perfoem the halo communications
      //     (1) use Send / Recv
      //     (2) use Isend / Irecv
      //         --> can you overlap communication and compution in this way?
      


      /* --------------------------------------  */

      /* --------------------------------------
      * HALO EXCHANGE START (Send/Recv)  --- BLOCKING VERSION ---
      * ---------------------------------------  */

      /*
      MPI_Status status;
      double t0, comm_time_north_send=0, comm_time_north_recv=0, comm_time_south_send=0, comm_time_south_recv=0;
      double comm_time_east_send=0, comm_time_east_recv=0, comm_time_west_send=0, comm_time_west_recv=0;
      // NORTH
      if (neighbours[NORTH] != MPI_PROC_NULL) {
          printf("Rank %d: About to recv from NORTH\n", Rank); fflush(stdout);
          t0 = MPI_Wtime();
          MPI_Recv(&current_plane[0 * full_sizex + 1], sizex, MPI_DOUBLE, neighbours[NORTH], 0, myCOMM_WORLD, &status);
          comm_time_north_recv += MPI_Wtime() - t0;
          printf("Rank %d: Recv from NORTH\n", Rank); fflush(stdout);

          printf("Rank %d: About to send to NORTH\n", Rank); fflush(stdout);
          t0 = MPI_Wtime();
          MPI_Send(&current_plane[1 * full_sizex + 1], sizex, MPI_DOUBLE, neighbours[NORTH], 0, myCOMM_WORLD);
          comm_time_north_send += MPI_Wtime() - t0;
          printf("Rank %d: Sent to NORTH\n", Rank); fflush(stdout);
      }
      // SOUTH
      if (neighbours[SOUTH] != MPI_PROC_NULL) {
          printf("Rank %d: About to recv from SOUTH\n", Rank); fflush(stdout);
          t0 = MPI_Wtime();
          MPI_Recv(&current_plane[(sizey + 1) * full_sizex + 1], sizex, MPI_DOUBLE, neighbours[SOUTH], 1, myCOMM_WORLD, &status);
          comm_time_south_recv += MPI_Wtime() - t0;
          printf("Rank %d: Recv from SOUTH\n", Rank); fflush(stdout);

          printf("Rank %d: About to send to SOUTH\n", Rank); fflush(stdout);
          t0 = MPI_Wtime();
          MPI_Send(&current_plane[sizey * full_sizex + 1], sizex, MPI_DOUBLE, neighbours[SOUTH], 1, myCOMM_WORLD);
          comm_time_south_send += MPI_Wtime() - t0;
          printf("Rank %d: Sent to SOUTH\n", Rank); fflush(stdout);
      }
      // EAST
      if (neighbours[EAST] != MPI_PROC_NULL) {
          printf("Rank %d: About to recv from EAST\n", Rank); fflush(stdout);
          t0 = MPI_Wtime();
          MPI_Recv(buffers[RECV][EAST], sizey, MPI_DOUBLE, neighbours[EAST], 2, myCOMM_WORLD, &status);
          comm_time_east_recv += MPI_Wtime() - t0;
          printf("Rank %d: Recv from EAST\n", Rank); fflush(stdout);

          printf("Rank %d: About to send to EAST\n", Rank); fflush(stdout);
          t0 = MPI_Wtime();
          MPI_Send(buffers[SEND][EAST], sizey, MPI_DOUBLE, neighbours[EAST], 2, myCOMM_WORLD);
          comm_time_east_send += MPI_Wtime() - t0;
          printf("Rank %d: Sent to EAST\n", Rank); fflush(stdout);
      }
      // WEST
      if (neighbours[WEST] != MPI_PROC_NULL) {
          printf("Rank %d: About to recv from WEST\n", Rank); fflush(stdout);
          t0 = MPI_Wtime();
          MPI_Recv(buffers[RECV][WEST], sizey, MPI_DOUBLE, neighbours[WEST], 3, myCOMM_WORLD, &status);
          comm_time_west_recv += MPI_Wtime() - t0;
          printf("Rank %d: Recv from WEST\n", Rank); fflush(stdout);

          printf("Rank %d: About to send to WEST\n", Rank); fflush(stdout);
          t0 = MPI_Wtime();
          MPI_Send(buffers[SEND][WEST], sizey, MPI_DOUBLE, neighbours[WEST], 3, myCOMM_WORLD);
          comm_time_west_send += MPI_Wtime() - t0;
          printf("Rank %d: Sent to WEST\n", Rank); fflush(stdout);
      }
      double comm_time_this_iter = comm_time_north_send + comm_time_north_recv +
                                   comm_time_south_send + comm_time_south_recv +
                                   comm_time_east_send + comm_time_east_recv +
                                   comm_time_west_send + comm_time_west_recv;
      total_comm_time += comm_time_this_iter;
      */

      /* ---------------------------------------------------- */
      /* HALO EXCHANGE START (Non-Blocking)      */
      /* ---------------------------------------------------- */


      // [B] Post all non-blocking communications
      // Post all RECEIVES first
      if (neighbours[NORTH] != MPI_PROC_NULL) {
          MPI_Irecv(&current_plane[0 * full_sizex + 1], sizex, MPI_DOUBLE, neighbours[NORTH], 1, myCOMM_WORLD, &requests[request_count]);
          request_count++;
      }
      if (neighbours[SOUTH] != MPI_PROC_NULL) {
          MPI_Irecv(&current_plane[(sizey + 1) * full_sizex + 1], sizex, MPI_DOUBLE, neighbours[SOUTH], 0, myCOMM_WORLD, &requests[request_count]);
          request_count++;
      }
      if (neighbours[EAST] != MPI_PROC_NULL) {
          MPI_Irecv(buffers[RECV][EAST], sizey, MPI_DOUBLE, neighbours[EAST], 3, myCOMM_WORLD, &requests[request_count]);
          request_count++;  
      }
      if (neighbours[WEST] != MPI_PROC_NULL) {
          MPI_Irecv(buffers[RECV][WEST], sizey, MPI_DOUBLE, neighbours[WEST], 2, myCOMM_WORLD, &requests[request_count]);
          request_count++;
      }

      // Then post all SENDS
      if (neighbours[NORTH] != MPI_PROC_NULL) {
          MPI_Isend(&current_plane[1 * full_sizex + 1], sizex, MPI_DOUBLE, neighbours[NORTH], 0, myCOMM_WORLD, &requests[request_count]);
          request_count++;
      }
      if (neighbours[SOUTH] != MPI_PROC_NULL) {
          MPI_Isend(&current_plane[sizey * full_sizex + 1], sizex, MPI_DOUBLE, neighbours[SOUTH], 1, myCOMM_WORLD, &requests[request_count]);
          request_count++;
      }
      if (neighbours[EAST] != MPI_PROC_NULL) {
          MPI_Isend(buffers[SEND][EAST], sizey, MPI_DOUBLE, neighbours[EAST], 2, myCOMM_WORLD, &requests[request_count]);
          request_count++;
      }
      if (neighbours[WEST] != MPI_PROC_NULL) {
          MPI_Isend(buffers[SEND][WEST], sizey, MPI_DOUBLE, neighbours[WEST], 3, myCOMM_WORLD, &requests[request_count]);
          request_count++;
      }
      // --- Wait for all communications to complete ---
      // MPI_Waitall will pause here until every single Isend and Irecv is finished.
      MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);


      // [C] Unpack the received East/West data
      if (neighbours[WEST] != MPI_PROC_NULL) {
          for (int i = 0; i < sizey; i++) current_plane[(i + 1) * full_sizex + 0] = buffers[RECV][WEST][i];
      }
      if (neighbours[EAST] != MPI_PROC_NULL) {
          for (int i = 0; i < sizey; i++) current_plane[(i + 1) * full_sizex + sizex + 1] = buffers[RECV][EAST][i];
      }

      /*
      * --------------------------------------
      * HALO EXCHANGE END
      * --------------------------------------
      */


      total_comm_time += (MPI_Wtime() - section_start_time);

      section_start_time = MPI_Wtime();

      /* update grid points */
      
      update_plane( periodic, N, &planes[current], &planes[!current] );

      total_comp_time += (MPI_Wtime() - section_start_time);

      /* output if needed */
      if ( output_energy_stat_perstep ) output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
      /* swap plane indexes for the new iteration */
      current = !current;
      
    }
  
  t_end = MPI_Wtime();
  
  if (Rank == 0) {
      printf("\n--- Timing Results ---\n");
      printf("Total loop time: %f seconds\n", t_end - t_start);
      printf("Time spent in communication: %f seconds\n", total_comm_time);
      printf("Time spent in computation: %f seconds\n", total_comp_time);
      printf("----------------------\n\n");

      // Calculate per-thread timing statistics
      double min_comp_time = 1e9, max_comp_time = 0.0, avg_comp_time = 0.0;
      double total_thread_time = 0.0;
      
      printf("\n--- Per-Thread Timing Analysis ---\n");
      printf("Number of OpenMP threads: %d\n", g_n_omp_threads);
      printf("Number of iterations: %d\n", Niterations);
      printf("\n");
      
      for (int i = 0; i < g_n_omp_threads; i++) {
            double thread_time = g_per_thread_comp_time[i];
            total_thread_time += thread_time;
            if (thread_time < min_comp_time) min_comp_time = thread_time;
            if (thread_time > max_comp_time) max_comp_time = thread_time;
            avg_comp_time += thread_time;
            printf("Thread %2d computation time: %f seconds\n", i, thread_time);
      }
      avg_comp_time /= g_n_omp_threads;
      
      printf("\n--- Thread Timing Statistics ---\n");
      printf("----------------------------------------------------------\n");
      printf("Min thread computation time: %f seconds\n", min_comp_time);
      printf("Max thread computation time: %f seconds\n", max_comp_time);
      printf("Avg thread computation time: %f seconds\n", avg_comp_time);
      printf("----------------------------------------------------------\n\n");
    }


  output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
  
  memory_release(planes, buffers);
  free(g_per_thread_comp_time);
  
  MPI_Finalize();
  return 0;
}


/* ==========================================================================
   =                                                                        =
   =   routines called within the integration loop                          =
   ========================================================================== */





/* ==========================================================================
   =                                                                        =
   =   initialization                                                       =
   ========================================================================== */


uint simple_factorization( uint, int *, uint ** );

int initialize_sources( int       ,
			int       ,
			MPI_Comm  *,
			uint      [2],
			int       ,
			int      *,
			vec2_t  ** );


int memory_allocate ( const int       *neighbours  ,
		      const vec2_t     N           ,
		            buffers_t *buffers_ptr ,
		            plane_t   *planes_ptr
		      )

{
    /*
      here you allocate the memory buffers that you need to
      (i)  hold the results of your computation
      (ii) communicate with your neighbours

      The memory layout that I propose to you is as follows:

      (i) --- calculations
      you need 2 memory regions: the "OLD" one that contains the
      results for the step (i-1)th, and the "NEW" one that will contain
      the updated results from the step ith.

      Then, the "NEW" will be treated as "OLD" and viceversa.

      These two memory regions are indexed by *plate_ptr:

      planew_ptr[0] ==> the "OLD" region
      plamew_ptr[1] ==> the "NEW" region


      (ii) --- communications

      you may need two buffers (one for sending and one for receiving)
      for each one of your neighnours, that are at most 4:
      north, south, east amd west.      

      To them you need to communicate at most mysizex or mysizey
      daouble data.

      These buffers are indexed by the buffer_ptr pointer so
      that

      (*buffers_ptr)[SEND][ {NORTH,...,WEST} ] = .. some memory regions
      (*buffers_ptr)[RECV][ {NORTH,...,WEST} ] = .. some memory regions
      
      --->> Of course you can change this layout as you prefer
      
     */

  if (planes_ptr == NULL ){
        // Print an error message to standard error.
        fprintf(stderr, "ERROR in memory_allocate: received an invalid NULL pointer for planes.\n");
        // Abort the entire MPI job, returning an error code of 1.
        MPI_Abort(MPI_COMM_WORLD, 1);
  }


  if (buffers_ptr == NULL ){
        // Print an error message to standard error.
        fprintf(stderr, "ERROR in memory_allocate: received an invalid NULL pointer for buffers.\n");
        // Abort the entire MPI job, returning an error code of 1.
        MPI_Abort(MPI_COMM_WORLD, 1);
   } 

  // ··················································
  // allocate memory for data

  const int sizex = planes_ptr[OLD].size[_x_];
  const int sizey = planes_ptr[OLD].size[_y_];
  
  // Check for valid domain sizes
  if (sizex <= 0 || sizey <= 0) {
      fprintf(stderr, "ERROR in memory_allocate: invalid domain size sizex=%d, sizey=%d\n", sizex, sizey);
      MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  // we allocate the space needed for the plane plus a contour frame
  // that will contains data form neighbouring MPI tasks
  unsigned int frame_size = (sizex + 2) * (sizey + 2);
  
  // Check for potential integer overflow
  if (frame_size == 0 || frame_size / (sizex + 2) != (sizey + 2)) {
      fprintf(stderr, "ERROR in memory_allocate: integer overflow in frame_size calculation\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
  }

  planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[OLD].data == NULL ){
        fprintf(stderr, "ERROR in memory_allocate: failed to allocate memory for OLD plane (size: %u doubles)\n", frame_size);
        MPI_Abort(MPI_COMM_WORLD, 1);
  }
  memset ( planes_ptr[OLD].data, 0, frame_size * sizeof(double) );

  planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[NEW].data == NULL ){
        fprintf(stderr, "ERROR in memory_allocate: failed to allocate memory for NEW plane (size: %u doubles)\n", frame_size);
        free(planes_ptr[OLD].data); // Clean up the already allocated memory
        MPI_Abort(MPI_COMM_WORLD, 1);
  }
  memset ( planes_ptr[NEW].data, 0, frame_size * sizeof(double) );

  // ··················································
  // buffers for north and south communication 
  // are not really needed
  //
  // in fact, they are already contiguous, just the
  // first and last line of every rank's plane
  //
  // you may just make some pointers pointing to the
  // correct positions
  //

  // or, if you preer, just go on and allocate buffers
  // also for north and south communications

  // ··················································
  // allocate buffers
  if (neighbours[EAST] != MPI_PROC_NULL) {
        buffers_ptr[SEND][EAST] = (double*)malloc(sizey * sizeof(double));
        if(buffers_ptr[SEND][EAST] == NULL) {
            fprintf(stderr, "ERROR in memory_allocate: failed to allocate EAST send buffer\n");
            free(planes_ptr[OLD].data);
            free(planes_ptr[NEW].data);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        buffers_ptr[RECV][EAST] = (double*)malloc(sizey * sizeof(double));
        if(buffers_ptr[RECV][EAST] == NULL) {
            fprintf(stderr, "ERROR in memory_allocate: failed to allocate EAST recv buffer\n");
            free(planes_ptr[OLD].data);
            free(planes_ptr[NEW].data);
            free(buffers_ptr[SEND][EAST]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
  }
  if (neighbours[WEST] != MPI_PROC_NULL) {
        buffers_ptr[SEND][WEST] = (double*)malloc(sizey * sizeof(double));
        if(buffers_ptr[SEND][WEST] == NULL) {
            fprintf(stderr, "ERROR in memory_allocate: failed to allocate WEST send buffer\n");
            free(planes_ptr[OLD].data);
            free(planes_ptr[NEW].data);
            if (neighbours[EAST] != MPI_PROC_NULL) {
                free(buffers_ptr[SEND][EAST]);
                free(buffers_ptr[RECV][EAST]);
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        buffers_ptr[RECV][WEST] = (double*)malloc(sizey * sizeof(double));
        if(buffers_ptr[RECV][WEST] == NULL) {
            fprintf(stderr, "ERROR in memory_allocate: failed to allocate WEST recv buffer\n");
            free(planes_ptr[OLD].data);
            free(planes_ptr[NEW].data);
            if (neighbours[EAST] != MPI_PROC_NULL) {
                free(buffers_ptr[SEND][EAST]);
                free(buffers_ptr[RECV][EAST]);
            }
            free(buffers_ptr[SEND][WEST]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
  }

  // ··················································

  
  return 0;
}



int memory_release ( plane_t *planes, buffers_t *buffers_ptr )
  
{
  if ( planes != NULL )
    {
      if ( planes[OLD].data != NULL )
	free (planes[OLD].data);
      
      if ( planes[NEW].data != NULL )
	free (planes[NEW].data);
    }

  // Safely free buffer memory - only free if pointers are not NULL
  if (buffers_ptr != NULL) {
      if (buffers_ptr[SEND][EAST] != NULL) {
          free(buffers_ptr[SEND][EAST]);
          buffers_ptr[SEND][EAST] = NULL;
      }
      if (buffers_ptr[RECV][EAST] != NULL) {
          free(buffers_ptr[RECV][EAST]);
          buffers_ptr[RECV][EAST] = NULL;
      }
      if (buffers_ptr[SEND][WEST] != NULL) {
          free(buffers_ptr[SEND][WEST]);
          buffers_ptr[SEND][WEST] = NULL;
      }
      if (buffers_ptr[RECV][WEST] != NULL) {
          free(buffers_ptr[RECV][WEST]);
          buffers_ptr[RECV][WEST] = NULL;
      }
  }
      
  return 0;
}



int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm )
{

  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy ( plane, &system_energy );
  
  MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
  
  if ( Me == 0 )
    {
      if ( step >= 0 )
	printf(" [ step %4d ] ", step ); fflush(stdout);

      
      printf( "total injected energy is %g, "
	      "system energy is %g "
	      "( in avg %g per grid point)\n",
	      budget,
	      tot_system_energy,
	      tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );
    }
  
  return 0;
}

int initialize ( MPI_Comm *Comm,
		 int      Me,                  // the rank of the calling process
		 int      Ntasks,              // the total number of MPI ranks
		 int      argc,                // the argc from command line
		 char   **argv,                // the argv from command line
		 vec2_t  *S,                   // the size of the plane
		 vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
		 int     *periodic,            // periodic-boundary tag
		 int     *output_energy_stat,
		 int     *neighbours,          // four-int array that gives back the neighbours of the calling task
		 int     *Niterations,         // how many iterations
		 int     *Nsources,            // how many heat sources
		 int     *Nsources_local,
		 vec2_t **Sources_local,
		 double  *energy_per_source,   // how much heat per source
		 plane_t *planes,
		 buffers_t *buffers
		 )
{
  int halt = 0;
  int ret;
  int verbose = 0;
  
  // ··································································
  // set deffault values

  (*S)[_x_]         = 10000;
  (*S)[_y_]         = 10000;
  *periodic         = 0;
  *Nsources         = 4;
  *Nsources_local   = 0;
  *Sources_local    = NULL;
  *Niterations      = 1000;
  *energy_per_source = 1.0;

  if ( planes == NULL ) {
    fprintf(stderr, "ERROR: planes pointer is NULL in initialize function\n");
    return 1;
  }

  // Initialize plane sizes to 0
  planes[OLD].size[_x_] = 0;
  planes[OLD].size[_y_] = 0;
  planes[NEW].size[_x_] = 0;
  planes[NEW].size[_y_] = 0;
  
  for ( int i = 0; i < 4; i++ )
    neighbours[i] = MPI_PROC_NULL;

  for ( int b = 0; b < 2; b++ )
    for ( int d = 0; d < 4; d++ )
      buffers[b][d] = NULL;
  
  // ··································································
  // process the commadn line
  // 
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:")) != -1)
      {
	switch( opt )
	  {
	  case 'x': (*S)[_x_] = (uint)atoi(optarg);
	    break;

	  case 'y': (*S)[_y_] = (uint)atoi(optarg);
	    break;

	  case 'e': *Nsources = atoi(optarg);
	    break;

	  case 'E': *energy_per_source = atof(optarg);
	    break;

	  case 'n': *Niterations = atoi(optarg);
	    break;

	  case 'o': *output_energy_stat = (atoi(optarg) > 0);
	    break;

	  case 'p': *periodic = (atoi(optarg) > 0);
	    break;

	  case 'v': verbose = atoi(optarg);
	    break;

	  case 'h': {
	    if ( Me == 0 )
	      printf( "\nvalid options are ( values btw [] are the default values ):\n"
		      "-x    x size of the plate [10000]\n"
		      "-y    y size of the plate [10000]\n"
		      "-e    how many energy sources on the plate [4]\n"
		      "-E    how many energy sources on the plate [1.0]\n"
		      "-n    how many iterations [1000]\n"
		      "-p    whether periodic boundaries applies  [0 = false]\n\n"
		      );
	    halt = 1; }
	    break;
	    
	    
	  case ':': printf( "option -%c requires an argument\n", optopt);
	    break;
	    
	  case '?': printf(" -------- help unavailable ----------\n");
	    break;
	  }
      }

    if ( opt == -1 )
      break;
  }

  if ( halt )
    return 1;
  
  
  // ··································································
  /*
   * here we should check for all the parms being meaningful
   *
   */

  // ...

  
  // ··································································
  /*
   * find a suitable domain decomposition
   * very simple algorithm, you may want to
   * substitute it with a better one
   *
   * the plane Sx x Sy will be solved with a grid
   * of Nx x Ny MPI tasks
   */

  vec2_t Grid;
  double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
  int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) );

  
 /* if ( dimensions == 1 )
    {
      if ( (*S)[_x_] >= (*S)[_y_] )
	Grid[_x_] = Ntasks, Grid[_y_] = 1;
      else
	Grid[_x_] = 1, Grid[_y_] = Ntasks;
    }
  else
    {
      int   Nf;
      uint *factors;
      uint  first = 1;
      ret = simple_factorization( Ntasks, &Nf, &factors );
      
      for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ )
	first *= factors[i];

      if ( (*S)[_x_] > (*S)[_y_] )
	Grid[_x_] = Ntasks/first, Grid[_y_] = first;
      else
	Grid[_x_] = first, Grid[_y_] = Ntasks/first;
    }*/

// --- Choose the most square 2D process grid (Px, Py) ---
int Px = 1, Py = Ntasks;
int min_diff = Ntasks; // Start with a large difference
for (int i = 1; i <= sqrt(Ntasks); i++) {
    if (Ntasks % i == 0) {
        int j = Ntasks / i;
        if (abs(i - j) < min_diff) {
            Px = i;
            Py = j;
            min_diff = abs(i - j);
        }
    }
}
Grid[_x_] = Px;
Grid[_y_] = Py;
(*N)[_x_] = Grid[_x_];
(*N)[_y_] = Grid[_y_];

  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_]
  // ··································································
  // my cooridnates in the grid of processors
  //
  int X = Me % Grid[_x_];
  int Y = Me / Grid[_x_];

  // ··································································
  // find my neighbours
  //

  if ( Grid[_x_] > 1 )
    {  
      if ( *periodic ) {       
	neighbours[EAST]  = Y*Grid[_x_] + (Me + 1 ) % Grid[_x_];
	neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
      
      else {
	neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
	neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
    }

  if ( Grid[_y_] > 1 )
    {
      if ( *periodic ) {      
	neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks; }

      else {    
	neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
	neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); }
    }

  // ··································································
  // the size of my patch
  //

  /*
   * every MPI task determines the size sx x sy of its own domain
   * REMIND: the computational domain will be embedded into a frame
   *         that is (sx+2) x (sy+2)
   *         the outern frame will be used for halo communication or
   */
  
  vec2_t mysize;
  uint s = (*S)[_x_] / Grid[_x_];
  uint r = (*S)[_x_] % Grid[_x_];
  mysize[_x_] = s + (X < r);
  s = (*S)[_y_] / Grid[_y_];
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + (Y < r);

  planes[OLD].size[_x_] = mysize[_x_];
  planes[OLD].size[_y_] = mysize[_y_];
  planes[NEW].size[_x_] = mysize[_x_];
  planes[NEW].size[_y_] = mysize[_y_];
  

  if ( verbose > 0 )
    {
      if ( Me == 0 ) {
	printf("Tasks are decomposed in a grid %d x %d\n\n",
		 Grid[_x_], Grid[_y_] );
	fflush(stdout);
      }

      MPI_Barrier(*Comm);
      
      for ( int t = 0; t < Ntasks; t++ )
	{
	  if ( t == Me )
	    {
	      printf("Task %4d :: "
		     "\tgrid coordinates : %3d, %3d\n"
		     "\tneighbours: N %4d    E %4d    S %4d    W %4d\n",
		     Me, X, Y,
		     neighbours[NORTH], neighbours[EAST],
		     neighbours[SOUTH], neighbours[WEST] );
	      fflush(stdout);
	    }

	  MPI_Barrier(*Comm);
	}
      
    }

  
  // ··································································
  // allocae the needed memory
  //
  ret = memory_allocate(neighbours, *N, buffers, planes);  

  // ··································································
  // allocae the heat sources
  //
  ret = initialize_sources( Me, Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local );
  
  
  return 0;  
}


uint simple_factorization( uint A, int *Nfactors, uint **factors )
/*
 * rought factorization;
 * assumes that A is small, of the order of <~ 10^5 max,
 * since it represents the number of tasks
 #
 */
{
  int N = 0;
  int f = 2;
  uint _A_ = A;

  while ( f < A )
    {
      while( _A_ % f == 0 ) {
	N++;
	_A_ /= f; }

      f++;
    }

  *Nfactors = N;
  uint *_factors_ = (uint*)malloc( N * sizeof(uint) );

  N   = 0;
  f   = 2;
  _A_ = A;

  while ( f < A )
    {
      while( _A_ % f == 0 ) {
	_factors_[N++] = f;
	_A_ /= f; }
      f++;
    }

  *factors = _factors_;
  return 0;
}


int initialize_sources( int       Me,
			int       Ntasks,
			MPI_Comm *Comm,
			vec2_t    mysize,
			int       Nsources,
			int      *Nsources_local,
			vec2_t  **Sources )

{

  srand48(time(NULL) ^ Me);
  int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
  
  if ( Me == 0 )
    {
      for ( int i = 0; i < Nsources; i++ )
	tasks_with_sources[i] = (int)lrand48() % Ntasks;
    }
  
  MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );

  int nlocal = 0;
  for ( int i = 0; i < Nsources; i++ )
    nlocal += (tasks_with_sources[i] == Me);
  *Nsources_local = nlocal;
  
  if ( nlocal > 0 )
    {
      vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
      for ( int s = 0; s < nlocal; s++ )
	{
	  helper[s][_x_] = 1 + lrand48() % mysize[_x_];
	  helper[s][_y_] = 1 + lrand48() % mysize[_y_];
	}

      *Sources = helper;
    }
  
  free( tasks_with_sources );

  return 0;
}

