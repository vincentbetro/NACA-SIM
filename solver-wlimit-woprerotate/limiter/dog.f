c=================================== TIMLIM ===================================
c
c Flux Limiter from paper by Tim Barth
c
c=============================================================================
      subroutine TIMLIM(nnodes,nface,qnode,qmax,qmin,phi,ncolor,
     1                  ncount,gradx,grady,x,y,fptr,nsface,nvface,
     2                  nfface,isface,ivface,ifface)
c
      dimension qnode(4,nnodes),qmax(4,nnodes),qmin(4,nnodes)
      dimension phi(4,nnodes)
      dimension gradx(4,nnodes),grady(4,nnodes),x(nnodes),y(nnodes)
c
      integer fptr(4,nface),ncount(1)
      integer isface(1),ivface(1),ifface(1)
c
      common/info/title(20),xmach,alpha,Re,dt,tot,rmstol,
     1            ntt,ivisc,irest,iflux,ileast,ncyc
      common/runge/rkc(5),cfl1,cfl2,iflim,ipv,nstage,iramp,nitfo,jupdate
      common/history/rms(3000),clw(3000),cdw(3000),cmw(3000),xres(3000)
c
c First loop over the edges and find the maximum and minimum
c
      do 1000 i = 1,nnodes
        qmax(1,i) = qnode(1,i)
        qmax(2,i) = qnode(2,i)
        qmax(3,i) = qnode(3,i)
        qmax(4,i) = qnode(4,i)
c
        qmin(1,i) = qnode(1,i)
        qmin(2,i) = qnode(2,i)
        qmin(3,i) = qnode(3,i)
        qmin(4,i) = qnode(4,i)
c
        phi(1,i) = 1.0
        phi(2,i) = 1.0
        phi(3,i) = 1.0
        phi(4,i) = 1.0
c
 1000 continue
      if(iflim.eq.0.or.ntt.lt.nitfo)return
c
      nstart = 1
      do 1010 i = 1,ncolor
CDIR$ IVDEP
        do 1020 n = nstart,nstart + ncount(i) - 1
          node1 = fptr(1,n)
          node2 = fptr(2,n)
c
c First check node1
c
          qmax(1,node1) = max(qmax(1,node1),qnode(1,node2))
          qmax(2,node1) = max(qmax(2,node1),qnode(2,node2))
          qmax(3,node1) = max(qmax(3,node1),qnode(3,node2))
          qmax(4,node1) = max(qmax(4,node1),qnode(4,node2))
c
          qmin(1,node1) = min(qmin(1,node1),qnode(1,node2))
          qmin(2,node1) = min(qmin(2,node1),qnode(2,node2))
          qmin(3,node1) = min(qmin(3,node1),qnode(3,node2))
          qmin(4,node1) = min(qmin(4,node1),qnode(4,node2))
c
c Now for node2
c
          qmax(1,node2) = max(qmax(1,node2),qnode(1,node1))
          qmax(2,node2) = max(qmax(2,node2),qnode(2,node1))
          qmax(3,node2) = max(qmax(3,node2),qnode(3,node1))
          qmax(4,node2) = max(qmax(4,node2),qnode(4,node1))
c
          qmin(1,node2) = min(qmin(1,node2),qnode(1,node1))
          qmin(2,node2) = min(qmin(2,node2),qnode(2,node1))
          qmin(3,node2) = min(qmin(3,node2),qnode(3,node1))
          qmin(4,node2) = min(qmin(4,node2),qnode(4,node1))
c
 1020   continue
      nstart = nstart + ncount(i)
 1010 continue
c
c Now we have found the max and min of surrounding nodes
c so let do the extrapolation to the face and "limit"
c the gradient so no new extrema are produced
c
      nstart = 1
      do 1030 i = 1,ncolor
CDIR$ IVDEP
        do 1040 n = nstart,nstart + ncount(i) - 1
          node1 = fptr(1,n)
          node2 = fptr(2,n)
c
          xmean = .5*(x(node1) + x(node2))
          ymean = .5*(y(node1) + y(node2))
c
c First do node 1
c
          rx = xmean - x(node1)
          ry = ymean - y(node1)
c
          q1 = qnode(1,node1) + gradx(1,node1)*rx + grady(1,node1)*ry
          q2 = qnode(2,node1) + gradx(2,node1)*rx + grady(2,node1)*ry
          q3 = qnode(3,node1) + gradx(3,node1)*rx + grady(3,node1)*ry
          q4 = qnode(4,node1) + gradx(4,node1)*rx + grady(4,node1)*ry
c
c Now check to see if this needs limiting
c
          if(q1.gt.qnode(1,node1))then
           temp = (qmax(1,node1) - qnode(1,node1))/(q1 - qnode(1,node1))
           phit = min(1.0,temp)
           phi(1,node1) = min(phi(1,node1),phit)
          else if(q1.lt.qnode(1,node1))then
           temp = (qmin(1,node1) - qnode(1,node1))/(q1 - qnode(1,node1))
           phit = min(1.0,temp)
           phi(1,node1) = min(phi(1,node1),phit)
          end if
c
          if(q2.gt.qnode(2,node1))then
           temp = (qmax(2,node1) - qnode(2,node1))/(q2 - qnode(2,node1))
           phit = min(1.0,temp)
           phi(2,node1) = min(phi(2,node1),phit)
          else if(q2.lt.qnode(2,node1))then
           temp = (qmin(2,node1) - qnode(2,node1))/(q2 - qnode(2,node1))
           phit = min(1.0,temp)
           phi(2,node1) = min(phi(2,node1),phit)
          end if
c
          if(q3.gt.qnode(3,node1))then
           temp = (qmax(3,node1) - qnode(3,node1))/(q3 - qnode(3,node1))
           phit = min(1.0,temp)
           phi(3,node1) = min(phi(3,node1),phit)
          else if(q3.lt.qnode(3,node1))then
           temp = (qmin(3,node1) - qnode(3,node1))/(q3 - qnode(3,node1))
           phit = min(1.0,temp)
           phi(3,node1) = min(phi(3,node1),phit)
          end if
c
          if(q4.gt.qnode(4,node1))then
           temp = (qmax(4,node1) - qnode(4,node1))/(q4 - qnode(4,node1))
           phit = min(1.0,temp)
           phi(4,node1) = min(phi(4,node1),phit)
          else if(q4.lt.qnode(4,node1))then
           temp = (qmin(4,node1) - qnode(4,node1))/(q4 - qnode(4,node1))
           phit = min(1.0,temp)
           phi(4,node1) = min(phi(4,node1),phit)
          end if
c
c Now for node 2
c
          rx = xmean - x(node2)
          ry = ymean - y(node2)
c
          q1 = qnode(1,node2) + gradx(1,node2)*rx + grady(1,node2)*ry
          q2 = qnode(2,node2) + gradx(2,node2)*rx + grady(2,node2)*ry
          q3 = qnode(3,node2) + gradx(3,node2)*rx + grady(3,node2)*ry
          q4 = qnode(4,node2) + gradx(4,node2)*rx + grady(4,node2)*ry
c
c Now check to see if this need limiting
c
          if(q1.gt.qnode(1,node2))then
           temp = (qmax(1,node2) - qnode(1,node2))/(q1 - qnode(1,node2))
           phit = min(1.0,temp)
           phi(1,node2) = min(phi(1,node2),phit)
          else if(q1.lt.qnode(1,node2))then
           temp = (qmin(1,node2) - qnode(1,node2))/(q1 - qnode(1,node2))
           phit = min(1.0,temp)
           phi(1,node2) = min(phi(1,node2),phit)
          end if
c
          if(q2.gt.qnode(2,node2))then
           temp = (qmax(2,node2) - qnode(2,node2))/(q2 - qnode(2,node2))
           phit = min(1.0,temp)
           phi(2,node2) = min(phi(2,node2),phit)
          else if(q2.lt.qnode(2,node2))then
           temp = (qmin(2,node2) - qnode(2,node2))/(q2 - qnode(2,node2))
           phit = min(1.0,temp)
           phi(2,node2) = min(phi(2,node2),phit)
          end if
c
          if(q3.gt.qnode(3,node2))then
           temp = (qmax(3,node2) - qnode(3,node2))/(q3 - qnode(3,node2))
           phit = min(1.0,temp)
           phi(3,node2) = min(phi(3,node2),phit)
          else if(q3.lt.qnode(3,node2))then
           temp = (qmin(3,node2) - qnode(3,node2))/(q3 - qnode(3,node2))
           phit = min(1.0,temp)
           phi(3,node2) = min(phi(3,node2),phit)
          end if
c
          if(q4.gt.qnode(4,node2))then
           temp = (qmax(4,node2) - qnode(4,node2))/(q4 - qnode(4,node2))
           phit = min(1.0,temp)
           phi(4,node2) = min(phi(4,node2),phit)
          else if(q4.lt.qnode(4,node2))then
           temp = (qmin(4,node2) - qnode(4,node2))/(q4 - qnode(4,node2))
           phit = min(1.0,temp)
           phi(4,node2) = min(phi(4,node2),phit)
          end if
 1040   continue
      nstart = nstart + ncount(i)
 1030 continue
c
      return
      end
