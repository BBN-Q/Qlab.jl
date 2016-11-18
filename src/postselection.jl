using KernelDensity

"""
    initByMsmt(data, numSegments, numMeas, indMeas, selectMeas, selectSign, threshold, docond)

Condition data on measurement results on the same channel

numSegments = number of distinct segments in the sequence
numMeas = number of meauserements per segments
indMeas = list of indeces of the msmt's to compare for fidelity
selectMeas = indeces of msms't to postselect on. All the msm'ts in a
segment will be postselected on this one
selectSign = select the sign for postselection. 1 = keep greater than;
-1 = keep smaller than
threshold
docond. If 1, do postselection
"""
function initByMsmt(data,numSegments,numMeas,indMeas,selectMeas,selectSign,threshold,docond)

  ind0=indMeas[1,:]
  ind1=indMeas[2,:]
  numShots = Int64(length(data)/(numMeas*numSegments))
  data=real(data)'
  datamat = splitdata(data,numShots,numMeas*numSegments)
  bins = linspace(minimum(minimum(datamat)), maximum(maximum(datamat)),500)

  PDFvec = zeros(length(bins),numMeas*numSegments)
  PDFvec_con = zeros(length(bins),numMeas*numSegments)

  for kk=1:numMeas*numSegments
    PDFvec[:,kk] = kde(datamat[:,kk], bins).density
  end

  opt_thr = zeros(length(ind0),1)
  fidvec_un = zeros(length(ind0),1)
  fidvec_con = zeros(length(ind0),1)
  err0 = zeros(length(ind0),1)
  err1 = zeros(length(ind0),1)
  data_con = zeros(numSegments, numMeas-length(selectMeas))


  for kk=1:size(PDFvec,2)
    PDFvec[:,kk]=PDFvec[:,kk]/sum(PDFvec[:,kk])
    PDFvec[:,kk]=PDFvec[:,kk]/sum(PDFvec[:,kk])
  end

  for kk=1:length(ind0)
    fidvec_un[kk] = 1-0.5*(1-0.5*sum(abs(PDFvec[:,ind0[kk]]-PDFvec[:,ind1[kk]])))
    indmaximum = indmax(abs(cumsum(PDFvec[:,ind0[kk]])-cumsum(PDFvec[:,ind1[kk]])))
    @printf("Optimum unconditioned fid. for segm. %d and %d = %.3f\n", ind0[kk], ind1[kk], fidvec_un[kk])
    tempvec0 = cumsum(abs(PDFvec[:,ind0[kk]]))
    err0[kk] = tempvec0[indmaximum]
    @printf("Error for |0> for segm. %d and %d = %.3f\n", ind0[kk], ind1[kk], tempvec0[indmaximum])
    tempvec1 = 1-cumsum(abs(PDFvec[:,ind1[kk]]))
    err1[kk] = tempvec1[indmaximum]
    @printf("Error for |1> for segm. %d and %d = %.3f\n", ind0[kk], ind1[kk], tempvec1[indmaximum])
    opt_thr[kk] = bins[indmaximum]
    @printf("Optimum threshold for segments %d and %d = %.4f\n", ind0[kk], ind1[kk], opt_thr[kk])
  end

  if(docond>0)
    for mm=1:numSegments
      for kk=1:length(selectMeas)
        ind = selectMeas[kk]
        for jj=1:numShots
          if(mm>4)
          end
          if selectSign == 1 && datamat[jj,numMeas*(mm-1)+ind] < threshold
            datamat[jj,numMeas*(mm-1)+1:numMeas*(mm-1)+numMeas]=NaN
          elseif selectSign == -1 && datamat[jj,numMeas*(mm-1)+ind] > threshold
            datamat[jj,numMeas*(mm-1)+1:numMeas*(mm-1)+numMeas]=NaN
          end
        end
      end
    end
  end

  for jj=1:numSegments
    thismeas=1;
    for kk=1:numMeas
      tempvec = filter!(vec->~isnan(vec) ,datamat[:,numMeas*(jj-1)+kk])
      PDFvec_con[:,(jj-1)*numMeas+kk] = kde(tempvec, bins).density

      #hist(~isnan(datamat[:,numMeas*(jj-1)+kk]).*datamat[:,numMeas*(jj-1)+kk], bins)[2];
      if(size(find(selectMeas.==kk),1)==0)
        data_con[jj,thismeas] = nanmean(datamat[:,numMeas*(jj-1)+kk])
        thismeas=thismeas+1
      end
    end
  end

  #normalize
  for kk=1:size(PDFvec_con,2)
    PDFvec_con[:,kk]=PDFvec_con[:,kk]/sum(PDFvec_con[:,kk])
    PDFvec_con[:,kk]=PDFvec_con[:,kk]/sum(PDFvec_con[:,kk])
  end


  for kk=1:length(ind0)
    fidvec_con[kk] = 1-0.5*(1-0.5*sum(abs(PDFvec_con[:,ind0[kk]]-PDFvec_con[:,ind1[kk]])))
    @printf("Optimum conditioned fid. for segm. %d and %d = %.3f\n", ind0[kk], ind1[kk], fidvec_con[kk])
  end
  return bins, PDFvec, PDFvec_con,data_con,opt_thr,datamat, fidvec_un, err0, err1, fidvec_con
end

function initByMsmt_2D(data,Anum,numSegments,numMeas_A,indMeas,selectMeas,selectSign,docond,numCal;threshold=NaN)
  #data = data qubits
  #data = single ancilla qubit. Data to be postselected on
  #Anum = index of ancilla measurement
  #numSegments = number of distinct segments in the sequence
  #numMeas_A = number of ancilla measurements per segments per channel (including those to be conditioned on)
  #indMeas = list of indeces of the msmt's to compare for fidelity. 1st dim = segments; 2nd dim: A channel. Note: removed comparison between multiple segments
  #selectMeas = indeces of msms't to postselect on. All the msm'ts in a
  #segment will be postselected on this one
  #selectSign = select the sign for postselection. 1 = keep greater than;
  #-1 = keep smaller than
  #docond. If 1, do postselection
  #numCal = number of calibration sequences. Currently, it assumes that they only have 1
  #measurement per qubit.
  #threshold. If not specified, uses the optimum to maximize fidelity

  numMeas_D=numMeas_A
  #the number of data measurements can be different from numMeas_A. For
  #instance, there may be multiple conditioning rounds before data
  #tomography. However, there is currently the constraint that all channels,
  #even if in different cards, acquire the same number of shots.

  ind0=indMeas[1,:]
  ind1=indMeas[2,:]
  nbins = 500

  if ~isa(Anum, Array)
    Anum = [Anum]
  end
  if ~isa(threshold, Array)
    threshold = repmat([threshold],length(Anum), 1)
  end

  data_A = zeros(length(data[1]), length(Anum))
  for (k, Aind) in enumerate(Anum)
    data_A[:, k] = real(data[Aind])
  end
  num_data_meas = length(data)-length(Anum)
  data_D = Dict()
  datamat = Dict()
  var_con = Dict()
  data_con = Dict()
  ind=1
  for kk = 1:num_data_meas
    while ind in Anum
      ind = ind+1
    end
    data_D[kk] = data[ind]
    ind+=1
  end

  numShots_A = Int64(floor(length(data_A[:, 1])/(numMeas_A*(numSegments-numCal)+numCal))) #assume same number of meas for all As
  for kk = 1:num_data_meas
    datamat[kk] = splitdata(real(data_D[kk]),numShots_A,numMeas_D*(numSegments-numCal)+numCal)
  end

  datamat_A = zeros(numShots_A, numMeas_A*(numSegments-numCal)+numCal, length(Anum))
  bins = Array(LinSpace{Float64},length(Anum))
  PDFvec = zeros(nbins,numMeas_A*(numSegments-numCal)+numCal, length(Anum))

  for Aind = 1:length(Anum)

    datamat_A[:, :, Aind] = splitdata(data_A[:, Aind],numShots_A,numMeas_A*(numSegments-numCal)+numCal)
    bins[Aind] = linspace(minimum(minimum(datamat_A[:, :, Aind])), maximum(maximum(datamat_A[:, :, Aind])),nbins)

    #PDFvec_con = zeros(length(bins),numMeas*numSegments);

    for kk=1:numMeas_A*(numSegments-numCal)+numCal
      PDFvec[:,kk,Aind] = kde(datamat_A[:,kk,Aind], bins[Aind]).density
    end
  end

  opt_thr = zeros(size(ind0))
  fidvec_un = zeros(size(ind0))
  #fidvec_con = zeros(length(ind0),1);
  err0 = zeros(size(ind0))
  err1 = zeros(size(ind0))
  datatemp = zeros(numSegments, numMeas_D)
  vartemp = zeros(numSegments, numMeas_D)
  #renormalize
  for Aind=1:length(Anum)
    for kk=1:size(PDFvec,2)
      PDFvec[:,kk,Aind]=PDFvec[:,kk,Aind]/sum(PDFvec[:,kk,Aind])
      PDFvec[:,kk,Aind]=PDFvec[:,kk,Aind]/sum(PDFvec[:,kk,Aind])
    end
  end

  for Aind = 1:length(Anum)
    @printf("Ancilla n.%d\n", Aind)
    #for kk=1:size(ind0,1)
    fidvec_un[Aind] = 1-0.5*(1-0.5*sum(abs(PDFvec[:,ind0[Aind],Aind]-PDFvec[:,ind1[Aind],Aind])))
    indmaximum = indmax(abs(cumsum(PDFvec[:,ind0[Aind],Aind])-cumsum(PDFvec[:,ind1[Aind],Aind])))
    @printf("Optimum unconditioned fid. for segm. %d and %d = %.3f\n", ind0[Aind], ind1[Aind], fidvec_un[Aind])
    tempvec0 = cumsum(abs(PDFvec[:,ind0[Aind],Aind]))
    err0[Aind] = tempvec0[indmaximum]
    @printf("Error for |0> for segm. %d and %d = %.3f\n", ind0[Aind], ind1[Aind], tempvec0[indmaximum])
    tempvec1 = 1-cumsum(abs(PDFvec[:,ind1[Aind],Aind]))
    err1[Aind] = tempvec1[indmaximum]
    @printf("Error for |1> for segm. %d and %d = %.3f\n", ind0[Aind], ind1[Aind], tempvec1[indmaximum])
    opt_thr[Aind] = bins[Aind][indmaximum]
    @printf("Optimum threshold for segments %d and %d = %.4f\n", ind0[Aind], ind1[Aind], opt_thr[Aind])
  end

  if docond>0
    if any(isnan(threshold))
      threshold  = opt_thr #uses one of the fidelity-maximizing values if unspecified
    end
    @printf("Length =%d\n", length(data_D))
    dataslice = datamat[1]
    for nn = 1:length(data_D)
      dataslice = datamat[nn]
      for Aind = 1:length(Anum)
        for mm=1:numSegments-numCal
          for kk=1:length(selectMeas)
            ind = selectMeas[kk]

            for jj=1:Int(numShots_A*numMeas_A/numMeas_D)
              if selectSign[kk, Aind] == 1 && datamat_A[jj,numMeas_A*(mm-1)+ind,Aind] < threshold[Aind]
                dataslice[jj,numMeas_A*(mm-1)+1:numMeas_A*(mm-1)+numMeas_A]=NaN
              elseif selectSign[kk, Aind] == -1 && datamat_A[jj,numMeas_A*(mm-1)+ind,Aind] > threshold[Aind]
                dataslice[jj,numMeas_A*(mm-1)+1:numMeas_A*(mm-1)+numMeas_A]=NaN
              end
            end

          end
        end
        datamat[nn] = dataslice
        frac = sum(sum(~isnan(dataslice[:,1:end-numCal])))/(size(dataslice,1)*(size(dataslice,2)-numCal))
        @printf("Fraction kept A%d = %.2f\n", Aind, frac)
      end

    end
  end

  for nn=1:length(data_D)
    dataslice = datamat[nn]
    for jj=1:numSegments
      thismeas=1
      for kk=1:numMeas_D
        if jj<numSegments-numCal+1
          datatemp[jj, thismeas] = nanmean(dataslice[:,numMeas_D*(jj-1)+kk])
          vartemp[jj, thismeas] = nanvar(dataslice[:, numMeas_D*(jj-1)+kk])
        else
          datatemp[jj, thismeas] = nanmean(dataslice[:,numMeas_D*(numSegments-numCal) + (jj-(numSegments-numCal))])
          vartemp[jj, thismeas] = nanvar(dataslice[:,numMeas_D*(numSegments-numCal) + (jj-(numSegments-numCal))]) #same cal. points for all data measurements
        end
        thismeas=thismeas+1
      end
    end

    data_con[nn] = datatemp[:]
    var_con[nn] = vartemp[:]

  end

  return (bins,PDFvec,data_con,var_con,opt_thr,datamat, fidvec_un, err0, err1)
end

function splitdata(data, numseg, numsets)
  datamat = zeros(numseg, numsets);
  for kk=1:numsets
    datamat[:,kk] = data[kk:numsets:end];
  end
  return datamat
end
