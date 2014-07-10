#include <Rcpp.h>

#include <sstream>
#include <cstdio>

#include <zlib.h>
#include <fcntl.h>

using namespace Rcpp;
using namespace std;

//' Write the results of an association test to a file
//' @param outfilename Output file name
//' @param indexfilename Index file name for the output
//' @param recordno An unsigned integer to use as a label for pvblock.
//' @param pvblock A data frame. The return value of \link[diseaseSims]{makePVblock}.
//' @param append  Boolean: append to output files or not?
//' @param gzip Boolean: write compressed data or not?
//' @details If append == TRUE, then file output is handled by low-level 
//' POSIX file locking via the C header
//' <fcntl.h>.  In practice, most file systems support this, but if you get
//' garbled output, yours may not.  The index file is locked, then the main
//' output is written and outfilename is closed.  The index data are then written,
//' and file locks are then released.
//'
//' If append == FALSE, then no file locking is done, which means
//' different R processes should write to different output files,
//' otherwise you'll be overwriting your data!!!!
// [[Rcpp::export]]
void writePVblock( const char * outfilename,
		   const char * indexfilename,
		   const unsigned & recordno,
		   DataFrame pvblock,
		   const bool & append = true,
		   const bool & gzip = false)
{
  NumericVector pos = pvblock["pos"],
    esize = pvblock["esize"],
    mfcontrols = pvblock["mfcontrols"],
    mfcases = pvblock["mfcases"],
    popfreq = pvblock["popfreq"],
    score = pvblock["score"];

  ostringstream out,error;
  out << "pos esize mfcontrols mfcases popfreq score\n";

  for( int i = 0 ; i < pvblock.nrows() ; ++i )//pos.size() ; ++i )
    {
      out << pos[i] << ' '
	  << esize[i] << ' '
	  << mfcontrols[i] << ' '
	  << mfcases[i] << ' '
	  << popfreq[i] << ' '
	  << score[i] << '\n';
    }

  if ( append )
    {
      struct flock index_flock;
      index_flock.l_type = F_WRLCK;/*Write lock*/
      index_flock.l_whence = SEEK_SET;
      index_flock.l_start = 0;
      index_flock.l_len = 0;/*Lock whole file*/

      FILE * index_fh = fopen(indexfilename,"a");
      int index_fd = fileno(index_fh);
      if ( index_fd == -1 ) 
	{ 
	  error << "ERROR: could not open " << indexfilename << '\n';
	  Rcpp::stop(error.str());
	}
      if (fcntl(index_fd, F_SETLKW,&index_flock) == -1) 
	{
	  error << "ERROR: could not obtain lock on " << indexfilename << '\n';
	  Rcpp::stop(error.str());
	}

      if (gzip)
	{
	  gzFile ofile = gzopen(outfilename,"a");
	  if( ofile == NULL )
	    {
	      error << "writePVblock error: could not open "
		    << outfilename << " for writing in append mode";
	      Rcpp::stop(error.str());
	   } 
	  int written = gzwrite( ofile, out.str().c_str(),out.str().size() );
	  gzclose(ofile);

	  //write data to index file AFTER writing
	  fprintf(index_fh,"%ul %d %d\n",recordno,written,pvblock.nrows());
	}
      else
	{	  
	  FILE * ofile = fopen(outfilename,"a");
	  if( ofile == NULL )
	    {
	      error << "writePVblock error: could not open "
		    << outfilename << " for writing in append mode\n";
	      Rcpp::stop(error.str());
	    }
	  //write data to index file BEFORE writing
	  fprintf(index_fh,"%ul %ld %d\n",recordno,ftell(ofile),pvblock.nrows());

	  fprintf(ofile,"%s",out.str().c_str());
	  fclose(ofile);
	}

      //release locks & close index file
      index_flock.l_type = F_UNLCK;
      if (fcntl(index_fd, F_UNLCK,&index_flock) == -1) 
	{
	  error << "ERROR: could not release lock on " << indexfilename << '\n';
	  Rcpp::stop(error.str());
	}
      fflush( index_fh );
      fclose(index_fh);
    }
  else
    {
      if( gzip )
	{
	  gzFile ofile = gzopen(outfilename,"w");
	  if( ofile == NULL )
	    {
	      error << "writePVblock error: could not open "
		    << outfilename << " for writing in append mode";
	      Rcpp::stop(error.str());
	    }
	  int written = gzwrite( ofile, out.str().c_str(), out.str().size() );
	  gzclose(ofile);

	  FILE * idx = fopen(indexfilename,"w");
	  if( idx == NULL )
	    {
	      error << "writePVblock error: could not open "
		    << indexfilename << " for writing in append mode";
	      Rcpp::stop(error.str());
	    }
	  fprintf( idx,"%ul %d %d\n",recordno,written,pvblock.nrows());
	  fclose(idx);
	}
      else
	{
	  FILE * idx = fopen(indexfilename,"w");
	  if( idx == NULL )
	    {
	      error << "writePVblock error: could not open "
		    << indexfilename << " for writing in append mode";
	      Rcpp::stop(error.str());
	    }
	  fprintf( idx,"%ul 0 %d\n",recordno,pvblock.nrows() );
	  fclose(idx);

	  idx = fopen(outfilename,"w");
	  if( idx == NULL )
	    {
	      error << "writePVblock error: could not open "
		    << outfilename << " for writing in append mode";
	      Rcpp::stop(error.str());
	    }
	  fprintf(idx,"%s",out.str().c_str());
	  fclose(idx);
	}
    }
  return;
}
