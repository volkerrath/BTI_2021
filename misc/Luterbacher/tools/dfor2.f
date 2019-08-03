       program dfor
c
chelpini
c_______________________________________________________________
c
c     Programa que pasa de formato extra a ascii (6(1x,e12.5))
c   
c     Files:
c
c      open(10,file=in,form='unformatted')
c      open(11,file=out)
c  
c     I/O files:
c     I: extra file
c     O: ascii file
c
c________________________________________________________________
chelpfin
c
      parameter(nr=4)
      integer ia(nr),nx,ny
      parameter (nx=130,ny=70)
      real f(nx,ny)
      character*64 in,out
      
      write (*,*)'Inputfile?' 
      read(*,150)in
c      write(*,*)'Outputfile?'
c      read(*,150)out
      out='temp.dat'
150   format(a64)

      open(10,file=in,form='unformatted')
c
      open(11,file=out)
      open(12, file='list.dat')
c
 10      read(10,end=1990)(ia(j),j=1,4)
         read(10)((f(i,k),i=1,nx),k=1,ny)
c         write(11,5010)ia
         write(12,'(4i10)')ia
         write(11,'(130(1x,e12.5))')((f(i,k),i=1,nx),k=1,ny)
      goto 10

 
 1990 close(10)

      stop
      end


