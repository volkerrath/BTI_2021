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
      integer ia(nr)
      real f(300000)
      character*64 in,out
      
      write (*,*)'Inputfile?' 
      read(*,150)in
      write(*,*)'Outputfile?'
      read(*,150)out
150   format(a64)

      open(10,file=in,form='unformatted')
c
      open(11,file=out)
c
 10      read(10,end=1990)(ia(j),j=1,4)
         read(10)(f(i1),i1=1,ia(4))
c         write(11,5010)ia
         write(*,5010)ia
         write(11,5020)(f(i1),i1=1,ia(4))
      goto 10

 5010 format(4i10)
 5020 format(6(1x,e12.5))
 5030 format(6(1x,e12.7))

 1990 close(10)

      stop
      end


