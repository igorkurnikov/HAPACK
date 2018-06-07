c-------------------------------------------------------
      subroutine gtoc(mgrid,oldmid,scale,g,c)
c
c converts grid to real coordinates
c
      dimension g(3),c(3)
	dimension oldmid(3)
	goff = (mgrid + 1.)/2.
      do 9000 i = 1,3
         c(i) = (g(i) - goff)/scale + oldmid(i)
9000	continue
      return
      end

c-------------------------------------------------------
      subroutine ctog(mgrid,oldmid,scale,c,g)
c
c converts real to grid coordinates
c
c---------------------------------------------------
c
      dimension g(3),c(3)
	dimension oldmid(3)

	goff = (mgrid + 1.)/2.
      do 9000 i = 1,3
         g(i) = (c(i) - oldmid(i))*scale + goff
9000	continue
      return
      end
