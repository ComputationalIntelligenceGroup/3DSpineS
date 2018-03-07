#' 3D virtual spine
#'
#' Generate the 3D representation of the virtual spines
#'
#' @param newSpines a matrix where each row is a virtual spine
#' @param idx is number of the row number of the spine to representate
#'
#' @return rgl object 3D mesh3d
#'
#' @examples
#' newSpines <- spineSampling(model,nSpines=5,cluster=1,seed=1)
#' mesh <- simulation3Dmesh(newSpines,2,iterations=4)
#' shade3d(mesh,col="red")
#'
#' @export
simulation3Dmesh<-function(newSpines,idx,iterations = 4)
{
  ellipsesSpine<-simulation3Dellipses(newSpines,idx)
  spine<-simulation3Dfaces(ellipsesSpine$ellipses,ellipsesSpine$skeleton)

  mesh<-list()
  mesh$vb<-rbind(t(spine$vertices),1)
  mesh$it<-t(spine$faces)
  class(mesh)<-"mesh3d"
  subdivision<-vcgSubdivide(mesh,type="Loop",looptype="continuity",silent=TRUE,iterations = iterations)
  return(subdivision)
}

#' Simulate new spines
#'
#' Simulate virtual spines from the mixture of multivariate Gaussians
#'
#' @param model is a Mclust object
#' @param nSpines is an integer that denotes the number of spines that must be simulated
#' @param cluster is an integer value that denotes the cluster from it should be sampled new spines, when it is 0 means that it should be sampled from prior
#' @param seed is an integer value used to obtain reproducible results
#'
#' @return a dataframe with then nSpines simulated spines
#'
#' @examples
#' newSpines <- spineSampling(model,nSpines=100,cluster=0,seed=1)
#'
#' @export
spineSampling<-function(model, nSpines=100, cluster=0,seed=1)
{
  if(class(model)!="Mclust")
  {
    stop("Parameter model must be an Mclust object")
  }

  if(!is.numeric(nSpines) | nSpines < 1)
  {
    stop("Number of samples must be an positive integer")
  }

  if(!is.numeric(cluster) | cluster < 0 | cluster > model$G)
  {
    stop(paste("Parameter cluster must be a number between 0 and",model$G))
  }

  if(!is.numeric(seed))
  {
    stop("Parameter seed must be a numeric value")
  }
  simulatedSpines<-matrix(0, ncol=ncol(model$data), nrow=nSpines)
  numValidSpines<-0
  i=0

  h_idx<-grep("^h\\_*",colnames(model$data))
  ellipse_axes_idx<-grep("B\\_*",colnames(model$data))

  if(model$scaled)
  {
    suma<-attributes(model$data)[[3]]
    cociente<-attributes(model$data)[[4]]

    theta_idx<-grep("theta\\_*",colnames(model$data))
    PCA_azimuth_position <- grep("inst\\_Phi\\_*", colnames(model$data))
    PCA_theta_position <- grep("inst\\_Theta\\_", colnames(model$data))
    angular_var<-c(theta_idx,PCA_theta_position,PCA_azimuth_position)

    recover<-rep(1,ncol(model$data))
    recover[angular_var]<-1/100000
    data<-sweep(model$data,2,recover,"*")
    data<-sweep(sweep(data,2,cociente,"*"),2,suma,"+")
    data<-data.frame(data,cluster=model$classification)
  }else{
    data<-data.frame(model$data,cluster=model$classification)
  }

   while(numValidSpines<nSpines)
   {
     sampledSpines <- sampleSpines(model, nSpines, cluster, seed+i)
     j=1
     while(numValidSpines<nSpines & j<nrow(sampledSpines))
     {
       if(all(sampledSpines[j,c(h_idx,ellipse_axes_idx)]>0))
       {
         spineEllipses <- simulation3Dellipses(sampledSpines,j)
         if(validateSpine(append(spineEllipses$ellipses,list(spineEllipses$skeleton[nrow(spineEllipses$skeleton),])),data,cluster))
         {
           numValidSpines<-numValidSpines + 1
           simulatedSpines[numValidSpines,]<-sampledSpines[j,]
         }
       }
       j<-j+1
     }
       i<-i+1
   }
  colnames(simulatedSpines)<-colnames(model$data)

  return(simulatedSpines)
}

#' Validate morphology of spines
#'
#' Function that validates if a spine is morphologically feasible checking the overlaping of the ellipses
#'
#' @param ellipses is a list of ellipses that defines the skeleton of the spine
#'
#' @return bool is true if ellipse is correct or false if the spine is morphologically incorrect
#' @noRd
validateSpine<-function(ellipses,data,cluster)
{
  #Simulate spine is incorrect if the angle between two elipses is more than 60º
  PCA <- as.numeric(prcomp(ellipses[[2]])$rotation[,3])
  for(i in 3:(length(ellipses)-1))
  {
    PCA_before<-PCA
    PCA <- as.numeric(prcomp(ellipses[[i]])$rotation[,3])

    if(acos(abs(PCA%*%PCA_before))>(0.34))
    {
      return(FALSE)
    }
  }

  #It is also wrong if two elipses are secants. Check if all the points
  for(i in 3:(length(ellipses)-1))
  {
    PCA <- prcomp(ellipses[[i]])
    media<-colMeans(ellipses[[i]])

    ellipse<-translate3d(ellipses[[i-1]],-media[1],-media[2],-media[3])

    if(abs(sum(sign(ellipse%*%PCA$rotation[,3])))<360)
    {
      return(FALSE)
    }
  }

  data_cluster<-data[which(data$cluster==cluster),]

  min_th<-apply(data_cluster,2,min)[grep("^V\\_",colnames(model$data))]
  max_th<-apply(data_cluster,2,max)[grep("^V\\_",colnames(model$data))]

  sim_vol<-matrix(0,nrow=1,ncol=(length(ellipses)-1))

  for(i in 2:length(ellipses))
  {
    convex<-convhulln(rbind(ellipses[[i-1]],ellipses[[i]]),"FA")
    sim_vol[i-1]<-convex$vol
  }

  if(!all(data.table::between(sim_vol,min_th,max_th)))
  {
    return(FALSE)
  }

  mesh<-list()
  mesh$vb<-t(cbind(rbind(ellipses[[length(ellipses)-2]],ellipses[[length(ellipses)-1]]),1))
  mesh$it<-t(convhulln(t(mesh$vb)[,1:3],"FA")$hull)
  class(mesh)<-"mesh3d"

  x<-list()
  x$vb<-matrix(ellipses[[length(ellipses)]],ncol=1,nrow=3)
  x$normals<-matrix(c(0,0,1),nrow=3,ncol=1)
  class(x)<-"mesh3d"

  if(vcgRaySearch(x,mesh)$quality==1)
  {
    return(FALSE)
  }

  return(TRUE)
}

#' Sample new spines
#'
#' Given a model select a cluster and sample from it or, if no cluster is given, sample from the prior distribution of the mixture
#'
#' @param model a Mclust object
#' @param nSpines is an integer value that denotes the number of spines to simulate
#' @param cluster is an integer value that denotes the cluster from it should be sampled new spines, when it is 0 means that it should be sampled from prior
#' @param seed is an integer value used to obtain reproducible results
#'
#' @return a matrix nSpines X  64 that represented virtualized spines
#' @noRd
sampleSpines<-function(model, nSpines, cluster,seed){
  set.seed(seed)
  if(cluster>0){
    new_spines <- mvrnorm(nSpines,model$parameters$mean[,cluster], model$parameters$variance$sigma[,,cluster])
  }else{
    #Choose cluster randomly according to the prior probability
    priori_mass_function <- cumsum(model$parameters$pro)
    random_val <- runif(100,0,1)
    probs <- sweep(matrix(rep(random_val,model$G),ncol=model$G),2,priori_mass_function)

    lower<-apply(model$data,2,min)
    upper<-apply(model$data,2,max)

    #Generate new spines
    new_spines <- matrix(0,ncol=length(model$parameters$mean[,1]),nrow=nSpines)
    for(i in 1:nSpines)
    {
      #new_spines[i,]<-rtmvnorm(1,model$parameters$mean[,which(probs[i,]<0)[1]], model$parameters$variance$sigma[,,which(probs[i,]<0)[1]],lower,upper,algorithm="gibbs")
      new_spines[i,]<-mvrnorm(1, model$parameters$mean[,which(probs[i,]<0)[1]], model$parameters$variance$sigma[,,which(probs[i,]<0)[1]])
    }
  }


  curvature_idx<-grep("cos_phi_*",colnames(model$data))
  if(model$scaled)
  {
    #Recover original values before standardization
    suma<-attributes(model$data)[[3]]
    cociente<-attributes(model$data)[[4]]

    theta_idx<-grep("theta\\_*",colnames(model$data))
    PCA_azimuth_position <- grep("inst\\_Phi\\_*", colnames(model$data))
    PCA_theta_position <- grep("inst\\_Theta\\_", colnames(model$data))
    angular_var<-c(theta_idx,PCA_theta_position,PCA_azimuth_position)


    recover<-rep(1,ncol(new_spines))
    recover[angular_var]<-1/100000
    new_spines<-sweep(new_spines,2,recover,"*")
    new_spines<-sweep(sweep(new_spines,2,cociente,"*"),2,suma,"+")
  }else{
    data<-data.frame(model$data,cluster=model$classification)
  }

  #Correct cosAngle
  idx_cosAngle <- which(new_spines[,curvature_idx] < -1)
  new_spines[,curvature_idx][idx_cosAngle] <- new_spines[,curvature_idx][idx_cosAngle]+abs(new_spines[,curvature_idx][idx_cosAngle]+1)*2

  idx_cosAngle <- which(new_spines[,curvature_idx] > 1)
  new_spines[,curvature_idx][idx_cosAngle] <- new_spines[,curvature_idx][idx_cosAngle]-abs(new_spines[,curvature_idx][idx_cosAngle]-1)*2

  colnames(new_spines)<-colnames(model$data)

  return(new_spines)
}

#' Generate the 3D representation of the spine represented uniquely by the ellipses
#'
#' Generate the 3D representation of the spine represented uniquely by the ellipses
#'
#' @param newSpines a matrix where each entry represents the values of a simulated spines
#' @param idx an integer that determines which row of newSpines must me selected to simulate that spine
#'
#' @return a list of three elements, first one is a matrix with all the vertices of all the ellipses,
#' second is a list saving in each element a matrix with the vertices of a ellipse, third are the centroids of the ellipses.
#'
#' @noRd
simulation3Dellipses<-function(newSpines,idx)
{

  curve_points<-data.frame(x = 0, y = 0, z = 0)
  cartesianEllipses<-list()
  cartesianEllipses[[1]]<-c(0,0,0)

  #Column index of the variables
  height_position <- grep("^h\\_*", colnames(newSpines))
  theta_position <- grep("theta\\_*", colnames(newSpines))
  cos_phi_position <- grep("cos\\_phi\\_*", colnames(newSpines))
  B_r_position <- grep("B\\_r\\_*", colnames(newSpines))
  B_R_position <- grep("B\\_R\\_*",colnames(newSpines))
  PCA_azimuth_position <- grep("inst\\_Phi\\_*", colnames(newSpines))
  PCA_theta_position <- grep("inst\\_Theta\\_", colnames(newSpines))

  vector_length <- as.numeric(newSpines[idx, height_position])
  phi <- as.numeric(newSpines[idx, theta_position])
  cosAngle <- as.numeric(newSpines[idx, cos_phi_position])
  B_r <- as.numeric(newSpines[idx, B_r_position])
  B_R <- as.numeric(newSpines[idx, B_R_position])
  perp_vectors_sph <- matrix(as.numeric(newSpines[idx, c(PCA_azimuth_position, PCA_theta_position)]), ncol = 2)

  skeleton <- matrix(0, nrow = length(height_position), ncol = 3)
  skeleton[1,]=c(0,0,1)*vector_length[1];

  for(i in 1:length(cosAngle))
  {
    if(i==1)
    {
      x<-c(0,0,-1);
    }else{
      x<-skeleton[i-1,]-skeleton[i,];
      x<-x/normv(x);
    }

    y<-c(0,0,0);

    u <- vcrossp(x,c(1,0,0));
    u <- as.numeric(u) / normv(u);
    theta <- acos(x %*% c(1,0,0));
    rotationMatrix<-rotation_matrix(u,theta);

    xyPosition<-c(cosAngle[i],sqrt(1-cosAngle[i]^2),0);
    rotationMatrix2<-rotation_matrix(c(1,0,0),phi[i]);
    original<-t(rotationMatrix2%*%xyPosition);
    lastArray<-as.numeric(t(t(rotationMatrix)%*%t(original))) * vector_length[i+1];
    skeleton[i+1,]<-lastArray+skeleton[i,];
  }

  major_axis <- B_R
  minor_axis <- B_r

  perp_PCA_matrix <- spherical2Cartesian(perp_vectors_sph[ ,1], perp_vectors_sph[ ,2], 1)
  for (i in 1:nrow(perp_PCA_matrix))
  {
    u <- vcrossp(perp_PCA_matrix[i,],c(0,0,1));
    u <- as.numeric(u) / normv(u);
    rotation_angle <- acos(perp_PCA_matrix[i,] %*% c(0,0,1));
    rotationMatrix<-rotation_matrix(u,rotation_angle);

    centers<-as.numeric(rotationMatrix%*%skeleton[i,]);

    fit <- list()
    fit$angle <- 1
    fit$maj <- max(major_axis[i],minor_axis[i])
    fit$min <- min(major_axis[i],minor_axis[i])
    fit$center <- centers[1:2]

    points_ellipse <- get.ellipse(fit)
    x <- points_ellipse[,1]
    y <- points_ellipse[,2]
    z <- centers[3]

    u = vcrossp(c(cos(1), sin(1), 0),c(1,0,0));
    u = as.numeric(u) / normv(u);
    rotationZ<-rotation_matrix(u,1);

    rotated_ellipse<-matrix(c(x, y, matrix(0,nrow=length(x),ncol=1)),nrow=length(x),ncol=3);
    rotated_ellipse<-t(rotationZ%*%t(rotated_ellipse));

    rotated_ellipse[,3]<-centers[3];
    elipse<-t(t(rotationZ)%*%t(rotated_ellipse));

    recovered_ellipse <- t(t(rotationMatrix) %*% t(elipse));
    cartesianEllipses[[i+1]]<-recovered_ellipse
    curve_points <- rbind(curve_points, data.frame(x = recovered_ellipse[,1], y = recovered_ellipse[,2], z = recovered_ellipse[,3]))
  }

  curve_points <- rbind(curve_points,skeleton[nrow(skeleton),])
  return(list(ellipses=cartesianEllipses, vertices=curve_points,skeleton=skeleton))
}

#' Generate the surface between ellipses
#'
#' Generate the surface between each pair of consecutive ellipses
#'
#' @param ellipses is a list where each element is a matrix of dimension nx3 representing the X,Y,Z coordinates of the ellipse
#' @param skeleton is a matrix with the centroids of each ellipse
#'
#' @return a list of two elements where each element is a matrix, first element is the set of vertices defining the spine and the second the faces joining those vertices
#'
#' @noRd
simulation3Dfaces<-function(ellipses,skeleton)
{
  num_points<-nrow(ellipses[[2]])
  points<-rbind(ellipses[[1]],ellipses[[2]])
  faces<-cbind(1,seq(2,(nrow(points)-1),by=1),seq(3,nrow(points),by=1))
  faces<-rbind(faces,c(1,nrow(points),2))

  for(i in 2:(length(ellipses)-1))
  {
    first_point_1<-nrow(points)-num_points+1
    first_point_2<-nrow(points)+1

    idx_first_ellipse<-seq(first_point_1,first_point_1+num_points-1,by=1)
    idx_snd_ellipse<-seq(first_point_2,first_point_2+num_points-1,by=1)

    first_triang<-matrix(c(idx_first_ellipse,idx_snd_ellipse,idx_snd_ellipse+1),ncol=3,byrow=F)
    first_triang[nrow(first_triang),3]<-idx_snd_ellipse[1]

    snd_triang<-matrix(c(idx_snd_ellipse,idx_first_ellipse,idx_first_ellipse-1),ncol=3,byrow=F)
    snd_triang[1,3]<-idx_first_ellipse[length(idx_first_ellipse)]

    faces2<-rbind(first_triang,snd_triang)

    faces<-rbind(faces,first_triang,snd_triang)

    points<-rbind(points,ellipses[[i+1]])
  }

  first_point_1 <- nrow(points)-num_points+1
  first_point_2 <- nrow(points)+1

  temp_faces<-cbind(first_point_2,seq(first_point_1,nrow(points),by=1)+1,seq(first_point_1,nrow(points),by=1))
  temp_faces[nrow(temp_faces),2]<-first_point_1[1]

  faces<-rbind(faces,temp_faces)

  points<-rbind(points,skeleton[nrow(skeleton),])

  return(list(vertices = points, faces = faces))
}
