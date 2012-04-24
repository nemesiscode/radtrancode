      integer function bboff(npath,nlayin,nbin,ipath,ilay,ibin,
     1 	maxbbbin)
C     ***************************************************************
C     function to get correct position in new BBBIN array to make it
C     a more memory-efficient array,
C
C     Pat Irwin		27/9/99
C
C     ***************************************************************
      integer npath,nlayin(npath),nbin,ipath,ilay,ibin
      integer ioff,maxbbbin

      ioff = 0
      do 10 i1 = 1,ipath-1 
       ioff = ioff + nlayin(i1)*nbin
10    continue

      ioff = ioff + nbin*(ilay-1)+ibin
   
      if(ioff.gt.maxbbbin)then
       print*,'Error in bboff.f ioff > maxbbbin'
       print*,'ioff, maxbbbin = ',ioff,maxbbbin
       stop
      endif

      bboff = ioff

      return

      end
