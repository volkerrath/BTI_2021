       parameter (nz=301)
       double precision t(nz) ,p(nz),z(nz),r(nz),k(nz),c(nz)
       double precision pmpa,tk
       
       do i=1,nz
          z(i)=(i-1)*10.d0
          t(i)=z(i)*0.03d0
          p(i)=1.d5+z(i)*998.d0*9.81d0
	  tk=t(i)+273.15d0
	  pmpa=p(i)/1.d6
c	  write (*,*) pmpa,tk
	  r(i)=dbpt(pmpa,tk)
	  k(i)=tcobpt(pmpa,tk)/1000.d0
	  c(i)=cpbpt(pmpa,tk)*1000.d0
	  write (*,*)  z(i),t(i),p(i),r(i),k(i),c(i)
       end do
       end
