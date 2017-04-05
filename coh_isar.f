c     ---------------------------------------------------------------
c     programme coh_isar.f
c     moyenne des fichiers "pha" InSAR
c     v2
c     5/4/2017 PB
c
c     warning: 'advance' statement in two 'write' instruction is not understood
c       by all f77 compilers, but ok with gfortran
c     ---------------------------------------------------------------

      parameter (nlig=5000,ncol=5000)
      character fich_in*62,fich_out*62,cr
c     character ligne*90
      integer*1 ival,ilig(nlig,ncol)
      integer iflag1(nlig,ncol),iii,icn,ira,ifla,imoy,iscat,idel
      integer icok,icnok,tab(9,2),isca(nlig,ncol),irange(256),iva
      integer*8 igt

      write(*,*) 'coh_insar version 2 20170404'
      write(*,*) 'Attention: - Taille maxi de l''image 5000 x 5000'
      write(*,*) '           - Tous les fichiers "pha*.*" seront lus'
      write(*,*) '           - Ce programme cree deux fichiers :'
      write(*,*) '              - un fichier de coherence'
      write(*,*) '              - un fichier de zones OK pour InSAR'
      write(*,*) ' '
      write(*,111)
  111 format('Nombre de lignes et de colonnes : ', $)
      read(*,*) nnlig,nncol

      do 39 ii=1,nnlig
      do 39 jj=1,nncol
        iflag1(ii,jj)=0
   39   isca(ii,jj)=0

      do 130 kk=1,256
  130   irange(kk)=0

c     Cherche les fichiers de phase dans le repertoire courant
      call system("ls pha*.*raw > fich")

      write(*,110)
  110 format('Fichier de coherence a ecrire   : ', $)
      read(*,'(a62)') fich_out
      open(25,file='fich')

c     Boucle de recherche des pixels noirs (phase zero)
      icn=0
   52 read(25,'(a)',end=59) fich_in
      icn=icn+1
      open(21,file=fich_in,access='stream')
      cr=char(13)
      write(*,'(a8,i4,1x,a35,a)',advance='no') 'Fichier ',icn,fich_in,cr
      do 53 i12=1,nnlig
      do 53 j12=1,nncol
c       Range I1 : -128;127
        read(21) ival
c       Tunable - here the sea is at 129
        ira=ival+129
        irange(ira)=irange(ira)+1
        if(ira.eq.1) iflag1(i12,j12)=iflag1(i12,j12)+1
   53 continue
      close(21)
      goto 52
   59 write(*,*) 'Fin recherche pixels noirs             '

      open(35,file='histogram.txt')

      write(35,*) 'Histogramme global'
      itop=0
      do 140 kkk=1,256
        if(irange(kkk).ge.itop) then
        ikkk=kkk
        itop=irange(kkk)
        endif
  140 write(35,*) kkk,irange(kkk)
      write(*,*) 'Background pixel : ',ikkk-1

c     To be implemented later
c      gnuplot -e "plot 'histo.txt' every ::2 notitle; set term gif;
c     + set output 'histo.gif'; replot"

c     Pixels blacks four times or more are discarded
      icok=0
      icnn=3+icn/20
      write(*,*) 'Threshold of number of black pixels : ',icnn
      open(26,file='noinsar.raw',access='stream')
      do 54 i13=1,nnlig
      do 54 j13=1,nncol
        iii=iflag1(i13,j13)
c       Pixels black three times or less are considered OK
c       The flag is set at zero
        if(iii.le.icnn) then
        icok=icok+1
        iflag1(i13,j13)=0
        ival=0
        else
c       Plots land in the range 128-255 (e.g. -127-0)
c       Can be modified in case of large set of input data
        if(iii.ge.132) iii=131
        ival=iii+124
        endif
   54 write(26) ival

      write(*,*) 'Pixels OK    : ',icok
      icnok=nnlig*nncol-icok
      write(*,*) 'Pixels noirs : ',icnok
      close(26)
      close(25)
              
c     Open output file
      open (18,file=fich_out,access='stream')
      numfi=0
      open(5,file='fich')

c     Start of big loop
      igt=0
    2 continue
     
      read(5,'(a)',end=1) fich_in
      numfi=numfi+1
      open(17,file=fich_in,access='stream')

      do 30 i5=1,2
      do 30 j5=1,nncol
   30   read(17) ilig(i5,j5)

      do 10 i=3,nnlig

      do 40 j=1,nncol
   40   read(17) ilig(3,j)

      do 41 j6=2,nncol-1

      imoy=0
      iicnt=0
      do 62 i11=1,3
      do 62 j11=1,3
        ifla=iflag1(i-3+i11,j6-2+j11)
        ill=ilig(i11,j6-2+j11)
        tab((3*i11-3)+j11,1)=ill
        tab((3*i11-3)+j11,2)=ifla
        if(ifla.eq.0) then
        iicnt=iicnt+1
        imoy=imoy+ill
        endif
   62 continue
      if(iicnt.ne.0) imoy=imoy/iicnt

c     Correct the mean value (for cycle slips) 
      do 77 i44=1,9
        if(tab(i44,2).eq.0) then
        idel=tab(i44,1)-imoy
        if(idel.gt.128) tab(i44,1)=tab(i44,1)-256
        if(idel.lt.-128) tab(i44,1)=tab(i44,1)+256
        endif
   77 continue
 
      imoy=0
      iicnt=0
      do 78 i45=1,9
        if(tab(i45,2).eq.0) then
        imoy=imoy+tab(i45,1)
        iicnt=iicnt+1
        endif
   78 continue
      if(iicnt.ne.0) imoy=imoy/iicnt

      goto 180

c     Corrects the mean value a second time (for cycle slips again)
      do 177 i46=1,9
        if(tab(i46,2).eq.0) then
          idel=tab(i46,1)-imoy
          if(idel.gt.128) then
            tab(i46,1)=tab(i46,1)-256
          endif
          if(idel.lt.-128) then
            tab(i46,1)=tab(i46,1)+256
          endif
        endif
  177 continue

      imoy=0
      iicnt=0
      do 178 i47=1,9
        if(tab(i47,2).eq.0) then
          imoy=imoy+tab(i47,1)
          iicnt=iicnt+1
        endif
  178 continue
      if(iicnt.ne.0) imoy=imoy/iicnt

  180 continue

c     Calcul de l'ecart type
c     Probleme de l'ambiguite de phase a peu pres pris en compte
c     Pas necessaire de faire mieux
      iscat=0
      do 63 k=1,9
        if(tab(k,2).eq.0) then
          idel=tab(k,1)-imoy
          if(abs(idel).ge.128) then
            if(idel.ge.0) then
              idel=idel-256
            else
              idel=idel+256
            endif
          endif
          iscat=iscat+idel*idel
        endif
   63 continue

c     Donne une note 255 aux clusters de moins de 6 bons pixels
      if(iicnt.ge.6) then
        iscat=dint(sqrt(dble(iscat)/dble(iicnt)))
        igt=igt+iscat
      else
        iscat=255
      endif

c     Pour chaque pixel additionne les ecarts type des interferogrammes successifs
      isca(i-1,j6)=isca(i-1,j6)+iscat
   41 continue

c     Decale des deux lignes de travail avant lecture de la ligne suivante
      do 42 i7=1,2
      do 42 j7=1,nncol
   42   ilig(i7,j7)=ilig(i7+1,j7)

   10 continue

      write(*,'(a12,i3,1x,a35,a)',advance='no') 'File number ',numfi
     +,fich_in,cr
      close(17)
      goto 2

    1 close(5)

      igt=(igt/(nncol-2))/(nnlig-2)

c Grand total
      write(*,*) ' '
      write(*,*) 'igt = ',igt

      do 190 kk=1,256
  190 irange(kk)=0

c     Ecriture du fichier de coherence
      ivam=0
      do 4 i1=1,nnlig
      do 4 j1=1,nncol

c     L'ecart type dans le fichier est le double de l'ecart type moyen
c     De facon a donner la maximum de dynamique au fichier
c     Parametre multiplicatif reglable entre 2 et 4
c     Ecrete aussi les pixels les plus coherents
        iva=dint(dble(4*isca(i1,j1))/dble(numfi))-32
        if(iva.lt.0) iva=0
        if(iva.gt.255) iva=-1
        if(iva.gt.ivam) ivam=iva
        ival=iva
        imax=255

c     Elimine les points non valides
        if(iflag1(i1,j1).ne.0) then
          ival=imax
        endif
        if(isca(i1,j1).eq.0) then
          ival=imax
          iva=0
        endif
        if(ival.gt.255) then
          ival=imax
          iva=ival
        endif

        if(iflag1(i1,j1).eq.0) then
          irange(iva+1)=irange(iva+1)+1
        endif

    4 write(18) ival
      close(18)
      write(*,*) cr,'Range of coherence : ',ivam,'               '
      
      open(135,file='histogram_coherence.txt')
      write(135,*) 'Histogram coherence'
      itop=0
      do 240 kk1=1,256
        write(135,*) kk1,irange(kk1)
  240 continue

      write(*,*) 'Done'
      call system("rm fich")
      stop
      end
