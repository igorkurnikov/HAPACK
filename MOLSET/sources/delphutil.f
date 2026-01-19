c  Fortran utility subroutines for DELPHI
c 
      subroutine openpm(nx,ny,nz, xstart,xend, ystart,yend, zstart,zend,
     $                  fname )
      
      integer nx, ny, nz
      real    xstart,xend, ystart,yend, zstat,zend 
      character*(*) fname

      real    extx, exty, extz
      real    xang, yang, zang 
      integer intx, inty, intz
      integer ivary, nbyte, intdat
      
      character*60 toplbl
      
      nx = 0
      ny = 0
      nz = 0
      xstart = 0.0
      xend   = 0.0
      ystart = 0.0
      yend   = 0.0
      zstart = 0.0
      zend   = 0.0

      open( unit=14, err=300, file=fname, status= "old",
     $      form="unformatted")
      
      read(14, err=310)toplbl
      read( 14, err=310) ivary,nbyte,intdat,extx,exty,extz,
     1  xang,yang,zang,xstart,xend,ystart,yend,zstart,
     1  zend,intx,inty,intz
      
      nx= intx + 1
      ny= inty + 1
      nz= intz + 1

      xstart = xstart * extx
      xend   = xend   * extx
      ystart = ystart * exty
      yend   = yend   * exty
      zstart = zstart * extz
      zend   = zend   * extz


      return
      
300      write(*,*) " Error in loadmp() "
      write(*,*) " error opening file ",fname 
      return
310      write(*,*) " Error in loadmp() "
      write(*,*) " reading potential map parameters from file ", fname
      close(14) 
      return

      end

      subroutine loadpm(nx,ny,nz, phimap )
c
c  Read Delphi Potential Map
c
      integer nx, ny, nz 
      real    phimap(*)
            
      do k = 1, nz
            do j = 1, ny 
              read(14,err=320)(phimap( nx*ny*(k - 1) + ny*(j -1) + i), 
     $                                 i =1, nx )  
            enddo
      enddo
      
      close(14)

      return
      
320      write(*,*) " Error in loadmp() "
      write(*,*) " reading potential map from file ", fname 
      close(14)
      return

      end
