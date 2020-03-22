      subroutine cosdir(cor,nn,i,j,k,cdir)
      include 'jumna_data.inc'
      dimension cor(nn,3),cdir(9)
      ax=cor(j,1)-cor(i,1)
      bx=cor(j,2)-cor(i,2)
      cx=cor(j,3)-cor(i,3)
      pm=sqrt(ax*ax+bx*bx+cx*cx)
      ax=ax/pm
      bx=bx/pm
      cx=cx/pm
      a=cor(k,1)-cor(j,1)
      b=cor(k,2)-cor(j,2)
      c=cor(k,3)-cor(j,3)
      az=bx*c-cx*b
      bz=cx*a-ax*c
      cz=ax*b-bx*a
      pm=sqrt(az*az+bz*bz+cz*cz)
      az=az/pm
      bz=bz/pm
      cz=cz/pm
      ay=bz*cx-cz*bx
      by=cz*ax-az*cx
      cy=az*bx-bz*ax
      pm=sqrt(ay*ay+by*by+cy*cy)
      cdir(1)=ax
      cdir(2)=ay/pm
      cdir(3)=az
      cdir(4)=bx
      cdir(5)=by/pm
      cdir(6)=bz
      cdir(7)=cx
      cdir(8)=cy/pm
      cdir(9)=cz
      return
      end
