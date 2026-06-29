// modified by M.K., Feb,2015

//////
// functions for file IO of SU(3) gauge fields
// by M.K. May,2009

inline double su3_unitalization( mdp_complex* p_su3 )
{
  double norm1 , norm2;
  mdp_complex proj;
  mdp_matrix su3( p_su3 , 3,3 );
//normalize first vector
//norm = <h1|h1>
    norm1 = sqrt( abs2(su3(0,0)) + abs2(su3(0,1)) + abs2(su3(0,2)) );
    for( int i=0 ; i<3 ; i++ ) su3(0,i) /= norm1;

//normalize second vector
//first orthogonalize h2 = h2 - <h1|h2>*h1
//projection = <h1|h2>
    proj = conj( su3(0,0) ) * su3(1,0)
	+ conj( su3(0,1) ) * su3(1,1)
	+ conj( su3(0,2) ) * su3(1,2);
    //if( abs(proj)>1.e-15 ) cout << "P";
    for( int i=0 ; i<3 ; i++ ) su3(1,i) -= proj * su3(0,i);

//then norm = <h2|h2>
    norm2 = sqrt( abs2(su3(1,0)) + abs2(su3(1,1)) + abs2(su3(1,2)) );
    for( int i=0 ; i<3 ; i++ ) su3(1,i) /= norm2;

//construct the vector product
    su3(2,0) = conj( su3(0,1)*su3(1,2) - su3(0,2)*su3(1,1) );
    su3(2,1) = conj( su3(0,2)*su3(1,0) - su3(0,0)*su3(1,2) );
    su3(2,2) = conj( su3(0,0)*su3(1,1) - su3(0,1)*su3(1,0) );
/*
  double norm3;
  norm3  = fabs( sqrt( abs2(su3(0,0)) + abs2(su3(1,0)) + abs2(su3(2,0)) ) -1.);
  norm3 += fabs( sqrt( abs2(su3(0,1)) + abs2(su3(1,1)) + abs2(su3(2,1)) ) -1.);
  norm3 += fabs( sqrt( abs2(su3(0,2)) + abs2(su3(1,2)) + abs2(su3(2,2)) ) -1.);
  if( norm3 > 1.e-14 ) cout << "error";
*/
    return ( fabs(norm1-1.)>fabs(norm2-1.) ) ? norm1 : norm2;
}

bool GaugeSave_SU3_mult( gauge_field& U , string filename ,
			 long max_buffer_size = 1024 )
{
  if( U.nc!=3 ) error( "This method is designed only for Nc=3\n" );

  const int nvol_s_1 = U.lattice().nx[1] * U.lattice().nx[2] * U.lattice().nx[3] -1 ;
  const int nlt = U.lattice().nx[0] / div0;

  const long field_components = U.size_per_site() * 2 / 3;
  const long psize = field_components * sizeof(mdp_complex);
  //mdp << nvol_s_1 << endl;
  //mdp << field_components << endl;
  //mdp << psize << endl;

  mpi.barrier();
  double mytime=mpi.time();
  mdp << "Saving file " << filename
      << " (buffer = " << max_buffer_size << " sites)" << endl; 
  mpi.barrier();

  //#pragma omp parallel for  remove: 2018 Jan 27 MK
  for( int t=0 ; t<U.lattice().nx[0] ; t++ ){
    int x[4];
    x[0] = t;
    // modify: 2018 Jan 27 MK
    const int processIO = 0;
    char dummy = 0;

    //char logfilename[100];
    //sprintf( logfilename , "file_log/%d_%d" , t , ME );
    //ofstream logfile( logfilename );

    if(ME==processIO) {
      long *buffer_size=new long[Nproc];
      long *buffer_ptr =new long[Nproc];
      mdp_array<mdp_complex,3> large_buffer(Nproc,max_buffer_size,field_components);
      mdp_complex *short_buffer=new mdp_complex[field_components];
      for( int process=0; process<Nproc ; process++ ) buffer_ptr[process]=0;

      char filename_num[100];
      sprintf( filename_num , "%s_%d" , filename.c_str() , t );
#pragma omp critical
      cout << "[t/proc/node= " << t << "/"
	   << processIO << "/" << omp_get_thread_num() << "]";
      fflush(stdout);

      FILE *fp=fopen( filename_num , "wb+" );
      if(fp==0) error( "Unable to open file" );

      for( x[1]=0 ; x[1]<U.lattice().nx[1] ; x[1]++ )
	for( x[2]=0 ; x[2]<U.lattice().nx[2] ; x[2]++ )
	  for( x[3]=0 ; x[3]<U.lattice().nx[3] ; x[3]++ ){

	    const int process = U.lattice().where( x , 4 , U.lattice().nx );
	  
	    if((process!=NOWHERE) && (process!=processIO)) {
	      if(buffer_ptr[process]==0) {
		//logfile << "G:" << x[1] << " " << x[2] << " " << x[3] << " : "
		//<< process << endl;
		mpi.get(buffer_size[process], process);
		mpi.get(&(large_buffer(process,0,0)), 
			buffer_size[process]*field_components, process);
		//mpi.put( dummy , process );// for blocking
	      }
	      for(int k=0 ; k<field_components ; k++ )
		short_buffer[k]=large_buffer(process,buffer_ptr[process],k);
	      buffer_ptr[process]++;
	      if(buffer_ptr[process]==buffer_size[process]) buffer_ptr[process]=0; 
	    }
	    if(process==processIO) {
	      for( int mu=0 ; mu<4 ; mu++ )
		for( int k=0 ; k<6 ; k++ )
		  short_buffer[mu*6+k] = U.mdp_complex_field::operator()( U.lattice().local(x) , mu*9+k );
	    }
	    if(process!=NOWHERE){
	      if( fwrite(short_buffer, psize,1, fp)!=1 )
		error("probably out of disk space");
	    }
	  }
      delete[] buffer_size;
      delete[] buffer_ptr;
      delete[] short_buffer;
      fclose(fp);
    } else {
      long buffer_size = 0 , idx;
      long *local_index = new long[max_buffer_size];
      mdp_array<mdp_complex,2> local_buffer(max_buffer_size,field_components);
      mdp_request request;

      for( x[1]=0 ; x[1]<U.lattice().nx[1] ; x[1]++ )
	for( x[2]=0 ; x[2]<U.lattice().nx[2] ; x[2]++ )
	  for( x[3]=0 ; x[3]<U.lattice().nx[3] ; x[3]++ ){
	    const int idx_s = x[3] + x[2] * U.lattice().nx[3]
	      + x[1] * U.lattice().nx[3] * U.lattice().nx[2];
	    const int process = U.lattice().where( x , 4 , U.lattice().nx );

	    if(process==ME) {
	      local_index[buffer_size] = U.lattice().local(x);
	      buffer_size++;
	    }
	    if((buffer_size==max_buffer_size) || 
	       ( ( idx_s == nvol_s_1 ) && ( buffer_size>0 ) ) ) {
	      for( idx=0 ; idx<buffer_size ; idx++)
		for( int mu=0 ; mu<4 ; mu++ )
		  for( int k=0 ; k<6 ; k++ )
		    local_buffer(idx,mu*6+k) = U.mdp_complex_field::operator()( local_index[idx] , mu*9+k );
	      
	      //logfile << "P:" << x[1] << " " << x[2] << " " << x[3] << " : "
	      //<< processIO << endl;

	      mpi.put(buffer_size, processIO, request);
	      mpi.wait(request);
	      mpi.put(&(local_buffer(0,0)), buffer_size*field_components, 
		      processIO, request);
	      mpi.wait(request);

	      //mpi.get( dummy , processIO );// for blocking
	      buffer_size=0;
	    }
	  }
      delete[] local_index;
    }
  }
  mpi.barrier();
  if(ME==0 && mdp_shutup==false) {
      printf("\n Saving time: %f (sec)\n", mpi.time()-mytime);
      fflush(stdout);
  }
  mpi.barrier();
  return true;
}


bool GaugeLoad_SU3_mult( gauge_field& U , string filename ,
			 long max_buffer_size = 1024 )
{
  if( U.nc!=3 ) error( "This method is designed only for Nc=3\n" );

  const int nvol_s_1 = U.lattice().nx[1] * U.lattice().nx[2] * U.lattice().nx[3] -1 ;
  const int nlt = U.lattice().nx[0] / div0;

  const long field_components = U.size_per_site() * 2 / 3;
  const long psize = field_components*sizeof(mdp_complex);

  mpi.barrier();

  double mytime=mpi.time();

  mdp << "Loading file " << filename
      << " (buffer = " << max_buffer_size << " sites)" << '\n'; 

  //#pragma omp parallel for remove: 2018 Jan 27 MK
  for( int t=0 ; t<U.lattice().nx[0] ; t++ ){
    int x[4];
    x[0] = t;
    // modify: 2018 Jan 27 MK
    const int processIO = 0;

    if(ME==processIO) {
      long *buffer_size=new long[Nproc];
      mdp_array<mdp_complex,3> large_buffer(Nproc,max_buffer_size,field_components);
      mdp_complex *short_buffer=new mdp_complex[field_components];
      int process;
      mdp_request request;

      for(process=0; process<Nproc; process++) buffer_size[process]=0;

      char filename_num[100];
      sprintf( filename_num , "%s_%d" , filename.c_str() , t );

#pragma omp critical
      cout << "[t/proc/node= " << t << "/"
	   << processIO << "/" << omp_get_thread_num() << "]";
      fflush(stdout);
      FILE *fp=fopen( filename_num , "rb" );
      if(fp==0) error( "Unable to open file" );

      for( x[1]=0 ; x[1]<U.lattice().nx[1] ; x[1]++ )
	for( x[2]=0 ; x[2]<U.lattice().nx[2] ; x[2]++ )
	  for( x[3]=0 ; x[3]<U.lattice().nx[3] ; x[3]++ ){

	    const int idx_s = x[3] + x[2] * U.lattice().nx[3]
	      + x[1] * U.lattice().nx[3] * U.lattice().nx[2];

	    process = U.lattice().where( x , 4 , U.lattice().nx );
	    if(process!=NOWHERE) {
	      if( fread(short_buffer, psize, 1, fp)!=1 ) 
		error("unexpected end of file");
	    }
	    if((process!=NOWHERE) && (process!=processIO)) {
	      for( int k=0 ; k<field_components ; k++ )
		large_buffer(process,buffer_size[process],k)=short_buffer[k];
	      buffer_size[process]++;
	      if(buffer_size[process]==max_buffer_size) {
		mpi.put(&(large_buffer(process,0,0)), 
			max_buffer_size*field_components, process, request); 
		mpi.wait(request);
		buffer_size[process]=0;
	      }
	      if(idx_s==nvol_s_1) 
		for(process=0; process<Nproc; process++)
		  if((process!=ME) && 
		     (buffer_size[process]!=max_buffer_size) &&
		     (buffer_size[process]>0)) {
		    mpi.put(&(large_buffer(process,0,0)), 
			    buffer_size[process]*field_components, 
			    process, request); 
		    mpi.wait(request);
		  }
	    }
	    if(process==processIO) {
	      for( int mu=0 ; mu<4 ; mu++ ){
		for( int k=0 ; k<6 ; k++ ){
		  U.mdp_complex_field::operator()( U.lattice().local(x) , mu*9+k ) = short_buffer[mu*6+k];
		}
		su3_unitalization( &U[ 36 * U.lattice().local(x) + mu*9 ] );
	      }
	    }
	  }
      delete[] buffer_size;
      delete[] short_buffer;
      fclose(fp);
    } else {
      int  process;
      long buffer_size=0, idx;
      long *local_index=new long[max_buffer_size];
      mdp_array<mdp_complex,2> local_buffer(max_buffer_size,field_components);

      for( x[1]=0 ; x[1]<U.lattice().nx[1] ; x[1]++ )
	for( x[2]=0 ; x[2]<U.lattice().nx[2] ; x[2]++ )
	  for( x[3]=0 ; x[3]<U.lattice().nx[3] ; x[3]++ ){

	    const int idx_s = x[3] + x[2] * U.lattice().nx[3]
	      + x[1] * U.lattice().nx[3] * U.lattice().nx[2];

	    process = U.lattice().where( x , 4 , U.lattice().nx );
	    if(process==ME) {
	      local_index[buffer_size]=U.lattice().local(x);
	      buffer_size++;
	    }
	    if((buffer_size==max_buffer_size) || 
	       ((idx_s==nvol_s_1) && (buffer_size>0))) {
	      mpi.get(&(local_buffer(0,0)), buffer_size*field_components, processIO);
	      for(idx=0; idx<buffer_size; idx++)
		for( int mu=0 ; mu<4 ; mu++ ){
		  for( int k=0 ; k<6 ; k++ ){
		    U.mdp_complex_field::operator()( local_index[idx] , mu*9+k ) = local_buffer(idx,mu*6+k);
		  }
		  su3_unitalization( &U[ 36 * local_index[idx] ] + mu*9 );
		}
	      buffer_size=0;
	    }
	  }
      delete[] local_index;
    }
  }
  U.update();
 
  mpi.barrier();
  if(ME==0 && mdp_shutup==false) {
    printf("\n Loading time: %f (sec)\n", mpi.time()-mytime);
    fflush(stdout);
  }
  mpi.barrier();
  return true;
}


bool GaugeSave_SU3( gauge_field& U , string filename ,
		    int processIO = 0 , 
		    long max_buffer_size = 1024, 
		    bool save_header = true,
		    long skip_bytes = 0 )
{
  if( U.nc!=3 ) error( "This method is designed only for Nc=3\n" );

  long header_size=0;
  const long field_components = U.size_per_site() * 2 / 3;
  const long psize = field_components * sizeof(mdp_complex);
  long idx_gl, nvol_gl=U.lattice().nvol_gl, k;
  double mytime=mpi.time();
  char dummy = 0;

  U.header.reset();
  if(ME==processIO) {
    long *buffer_size=new long[Nproc];
    long *buffer_ptr =new long[Nproc];
    mdp_array<mdp_complex,3> large_buffer(Nproc,max_buffer_size,field_components);
    mdp_complex *short_buffer=new mdp_complex[field_components];
    int process;
    for(process=0; process<Nproc; process++) buffer_ptr[process]=0;
    cout << "Saving file " << filename 
      << " from process " << processIO 
      << " (buffer = " << max_buffer_size << " sites)" << '\n'; 
    fflush(stdout);
    FILE *fp=fopen(filename.c_str(), "wb+");
    if(fp==0) error("Unable to open file");      

    U.header.set_time();

    if(save_header) {
      header_size=sizeof(mdp_field_file_header);
      if(fseek(fp, skip_bytes, SEEK_SET) || 
        fwrite(&U.header, header_size, 1, fp)!=1)
	error("Unable to write file header");
    }

    skip_bytes+=header_size;

    int x[4];
    for( x[0]=0 ; x[0]<U.lattice().nx[0] ; x[0]++ )
    for( x[1]=0 ; x[1]<U.lattice().nx[1] ; x[1]++ )
    for( x[2]=0 ; x[2]<U.lattice().nx[2] ; x[2]++ )
    for( x[3]=0 ; x[3]<U.lattice().nx[3] ; x[3]++ ){
      idx_gl= x[3]
	+ x[2] * U.lattice().nx[3] 
	+ x[1] * U.lattice().nx[3] * U.lattice().nx[2] 
	+ x[0] * U.lattice().nx[3] * U.lattice().nx[2] * U.lattice().nx[1];
      process = U.lattice().where( x , 4 , U.lattice().nx );

      if((process!=NOWHERE) && (process!=processIO)) {
	if(buffer_ptr[process]==0) {
	  mpi.get(buffer_size[process], process);
	  mpi.get(&(large_buffer(process,0,0)), 
		  buffer_size[process]*field_components, process);
          mpi.put( dummy , process );// for blocking
 	}
	for(k=0; k<field_components; k++)
	  short_buffer[k]=large_buffer(process,buffer_ptr[process],k);
	buffer_ptr[process]++;
	if(buffer_ptr[process]==buffer_size[process]) buffer_ptr[process]=0; 
      }
      if(process==processIO) {
	for( int mu=0 ; mu<4 ; mu++ )
	  for( k=0 ; k<6 ; k++ )
	    short_buffer[mu*6+k] = U.mdp_complex_field::operator()( U.lattice().local(x) , mu*9+k );
      }
      if(process!=NOWHERE) {
        if(//fseek(fp, idx_gl*psize+skip_bytes, SEEK_SET) ||
	   fwrite(short_buffer, psize,1, fp)!=1)
	    error("probably out of disk space");
      }
    }
    delete[] buffer_size;
    delete[] buffer_ptr;
    delete[] short_buffer;
    fclose(fp);
  } else {
    int  process;
    long buffer_size=0, idx, idx_gl;
    long *local_index=new long[max_buffer_size];
    mdp_array<mdp_complex,2> local_buffer(max_buffer_size,field_components);
    mdp_request request;
    int x[4];
    for( x[0]=0 ; x[0]<U.lattice().nx[0] ; x[0]++ )
    for( x[1]=0 ; x[1]<U.lattice().nx[1] ; x[1]++ )
    for( x[2]=0 ; x[2]<U.lattice().nx[2] ; x[2]++ )
    for( x[3]=0 ; x[3]<U.lattice().nx[3] ; x[3]++ ){
      idx_gl= x[3]
	+ x[2] * U.lattice().nx[3] 
	+ x[1] * U.lattice().nx[3] * U.lattice().nx[2] 
	+ x[0] * U.lattice().nx[3] * U.lattice().nx[2] * U.lattice().nx[1];
      process = U.lattice().where( x , 4 , U.lattice().nx );
      if(process==ME) {
	local_index[buffer_size]=U.lattice().local(x);
	buffer_size++;
      }
      if((buffer_size==max_buffer_size) || 
	 ((idx_gl==nvol_gl-1) && (buffer_size>0))) {
	for(idx=0; idx<buffer_size; idx++)
	  for( int mu=0 ; mu<4 ; mu++ )
	    for( k=0 ; k<6 ; k++ )
	      local_buffer(idx,mu*6+k) = U.mdp_complex_field::operator()( local_index[idx] , mu*9+k );

	mpi.put(buffer_size, processIO, request);
	mpi.wait(request);
	mpi.put(&(local_buffer(0,0)), buffer_size*field_components, 
		processIO, request);
	mpi.wait(request);
        mpi.get( dummy , processIO );// for blocking
 	buffer_size=0;
      }
    }
    delete[] local_index;
  }
  if(ME==0 && mdp_shutup==false) {
      printf("\n... Saving time: %f (sec)\n", mpi.time()-mytime);
      fflush(stdout);
  }
  return true;
}




bool GaugeLoad_SU3( gauge_field& U , string filename ,
                    int processIO = 0 ,
                    long max_buffer_size = 1024,
                    bool load_header = true,
                    long skip_bytes = 0 )
{
  if( U.nc!=3 ) error( "This method is designed only for Nc=3\n" );

  long header_size=0;
  const long field_components = U.size_per_site() * 2 / 3;
  long idx_gl, nvol_gl=U.lattice().nvol_gl, k;
  long psize=field_components*sizeof(mdp_complex);
  double mytime=mpi.time();
  bool reversed_header_endianess=false;
  if(ME==processIO) {
    long *buffer_size=new long[Nproc];
    mdp_array<mdp_complex,3> large_buffer(Nproc,max_buffer_size,field_components);
    mdp_complex *short_buffer=new mdp_complex[field_components];
    int process;
    mdp_request request;

    for(process=0; process<Nproc; process++) buffer_size[process]=0;
    cout << "Loading file " << filename 
	 << " from process " << processIO 
	 << " (buffer = " << max_buffer_size << " sites)" << '\n'; 
    fflush(stdout);
    FILE *fp=fopen(filename.c_str(), "rb");
    if(fp==0) error("Unable to open file");
    //cout << "opened" << endl;

    int i;
    if(load_header) {
      mdp_field_file_header tmp_header;
      header_size=sizeof(mdp_field_file_header);
      if(fseek(fp, skip_bytes, SEEK_SET) || 
	 fread(&tmp_header, header_size, 1, fp)!=1) {
	fprintf(stderr,"mdp_field.load(): Unable to load file header\n");
	return false;
      }
      
      reversed_header_endianess=switch_header_endianess(tmp_header);
      if(tmp_header.endianess!=U.header.endianess) 
	fprintf(stderr, "Unrecognized endianess... trying to read anyway\n");

      if(tmp_header.ndim!=U.header.ndim) {
	fprintf(stderr,"mdp_field.load(): wrong ndim\n");
	return false;
      }
      for(i=0; i<U.lattice().ndim; i++)  
	if(tmp_header.box[i]!=U.header.box[i]) {
	  fprintf(stderr,"mdp_file.load(): wrong lattice size\n");
	  return false;
	}
      if(tmp_header.bytes_per_site!=U.header.bytes_per_site) {
	fprintf(stderr, "mdp_file.load(): wrong type of field\n");
	return false;
      }
      if(tmp_header.sites!=U.header.sites) {
	fprintf(stderr,"mdp_field.load(): wrong number of sites\n");
	return false;
      }
      U.header=tmp_header;
    }
    
    skip_bytes+=header_size;
    
    int x[4];
    for( x[0]=0 ; x[0]<U.lattice().nx[0] ; x[0]++ )
    for( x[1]=0 ; x[1]<U.lattice().nx[1] ; x[1]++ )
    for( x[2]=0 ; x[2]<U.lattice().nx[2] ; x[2]++ )
    for( x[3]=0 ; x[3]<U.lattice().nx[3] ; x[3]++ ){
      idx_gl= x[3]
	+ x[2] * U.lattice().nx[3] 
	+ x[1] * U.lattice().nx[3] * U.lattice().nx[2] 
	+ x[0] * U.lattice().nx[3] * U.lattice().nx[2] * U.lattice().nx[1];

      process = U.lattice().where( x , 4 , U.lattice().nx );
      if(process!=NOWHERE) {
	if(//fseek(fp, idx_gl*psize+skip_bytes, SEEK_SET) ||
	   fread(short_buffer, psize, 1, fp)!=1) 
	  error("unexpected end of file");
      }
      if((process!=NOWHERE) && (process!=processIO)) {
	for(k=0; k<field_components; k++)
	  large_buffer(process,buffer_size[process],k)=short_buffer[k];
	buffer_size[process]++;
	if(buffer_size[process]==max_buffer_size) {
	  mpi.put(&(large_buffer(process,0,0)), 
		  max_buffer_size*field_components, process, request); 
	  mpi.wait(request);
	  buffer_size[process]=0;
	}
	if(idx_gl==nvol_gl-1) 
	  for(process=0; process<Nproc; process++)
	    if((process!=ME) && 
	       (buffer_size[process]!=max_buffer_size) &&
	       (buffer_size[process]>0)) {
	      mpi.put(&(large_buffer(process,0,0)), 
		      buffer_size[process]*field_components, 
		      process, request); 
	      mpi.wait(request);
	    }
      }
      if(process==processIO) {
	for( int mu=0 ; mu<4 ; mu++ ){
	  for( k=0 ; k<6 ; k++ )
	    U.mdp_complex_field::operator()( U.lattice().local(x) , mu*9+k ) = short_buffer[mu*6+k];
	  su3_unitalization( &U[ 36 * U.lattice().local(x) + mu*9 ] );
	}
      }
    }
    delete[] buffer_size;
    delete[] short_buffer;
    fclose(fp);
  } else {
    int  process;
    long buffer_size=0, idx;
    long *local_index=new long[max_buffer_size];
    mdp_array<mdp_complex,2> local_buffer(max_buffer_size,field_components);

    int x[4];
    for( x[0]=0 ; x[0]<U.lattice().nx[0] ; x[0]++ )
    for( x[1]=0 ; x[1]<U.lattice().nx[1] ; x[1]++ )
    for( x[2]=0 ; x[2]<U.lattice().nx[2] ; x[2]++ )
    for( x[3]=0 ; x[3]<U.lattice().nx[3] ; x[3]++ ){
      idx_gl= x[3]
	+ x[2] * U.lattice().nx[3] 
	+ x[1] * U.lattice().nx[3] * U.lattice().nx[2] 
	+ x[0] * U.lattice().nx[3] * U.lattice().nx[2] * U.lattice().nx[1];

      process = U.lattice().where( x , 4 , U.lattice().nx );
      if(process==ME) {
	local_index[buffer_size]=U.lattice().local(x);
	buffer_size++;
      }
      if((buffer_size==max_buffer_size) || 
	 ((idx_gl==nvol_gl-1) && (buffer_size>0))) {
	mpi.get(&(local_buffer(0,0)), buffer_size*field_components, processIO);
	for(idx=0; idx<buffer_size; idx++)
	  for( int mu=0 ; mu<4 ; mu++ ){
	    for( k=0 ; k<6 ; k++ )
	      U.mdp_complex_field::operator()( local_index[idx] , mu*9+k ) = local_buffer(idx,mu*6+k);
	  su3_unitalization( &U[ 36 * local_index[idx] ] + mu*9 );
	  }
	buffer_size=0;
      }
    }
    delete[] local_index;
  }

  U.update();
  mpi.broadcast(reversed_header_endianess,processIO);
 
  if(ME==0 && mdp_shutup==false) {
    printf("... Loading time: %f (sec)\n", mpi.time()-mytime);
    fflush(stdout);
  }
  return true;
}
