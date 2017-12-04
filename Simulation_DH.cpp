#include <Rcpp.h>
#include <iomanip>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
DataFrame getSampleDH(List parameter, 
                      Nullable<NumericVector> sampleSize = R_NilValue)
{
  
  int ss;
  
  int ssi = as<int> (parameter["ssi"]);
  int rp = as<int> (parameter["rp"]);
  
  if(sampleSize.isNotNull()) {
    ss = as<int> (sampleSize.get());
  } else {
    ss = ssi * rp;
  }
  
  // Disease Prevalence
  double din = as<double> (parameter["din"]);
  double dinP = as<double> (parameter["dinP"]);
  
  // var
  double heritability = as<double> (parameter["heritability"]);
  double eR2 = 1 - heritability;
  // sd
  double sdG = sqrt(heritability);
  double sdE = sqrt(eR2);
  
  // methylation
  double meanM = as<double> (as<List> (parameter["methylation"])["meanM"]);
  double varM = as<double> (as<List> (parameter["methylation"])["varM"]);
  double rsqE = as<double> (as<List> (parameter["methylation"])["rsqE"]);
  double scaleE = sqrt(varM * rsqE) / sqrt(1 - heritability);
  double varME = varM * (1 - rsqE);
  double corME = as<double> (as<List> (parameter["methylation"])["corME"]);
  //sd
  double sdME = sqrt(varME);
  
  // double Tcutoff = R::qnorm(1.0 - din, 0, 1, true, false);
  double TcutoffD = R::qnorm(1.0 - ((1 - dinP) * din), 0, 1, true, false);
  double TcutoffH = R::qnorm((1 - dinP) * (1 - din), 0, 1, true, false);
  
  
  bool twin = as<bool> (parameter["twin"]);
  if (!twin) {
    //** not twin **//
    
    // data.frame output
    NumericVector gValue(2*ss);
    NumericVector eValue(2*ss);
    NumericVector mValue(2*ss);
    NumericVector bValue(2*ss);
    NumericVector liability(2*ss);
    NumericVector DH(2*ss);
    
    // case and control index
    int nCase = 0;
    int nControl = 0;
    
    Rcout << "Getting sample..." << std::endl;
    
    // Rcout << "\r"
    //       << std::fixed << std::setprecision(2)
    //       << "Getting sample (Case[" << (double (nCase) / ss * 100) << "%] : "
    //       << nCase << " of " << ss
    //       << ", Control[" << (double (nControl) / ss * 100) << "%]: "
    //       << nControl << " of " << ss << ")";
    // // output counter
    // int sCounter = 0;
    // // print when counter is
    // const int sCounterTh = 10000;
    
    // case and control still require
    bool rCase = true;
    bool rControl = true;
    while (rCase | rControl) {
      
      double ilg = R::rnorm(0, sdG);
      double ile = R::rnorm(0, sdE);
      
      double iLiability = ilg + ile;
      
      bool iDH = 0;
      int iIdx = -1;
      if(iLiability > TcutoffD) {
        if(nCase < ss) {
          iDH = 1;
          iIdx = nCase;
          ++nCase;
        }
        else {
          rCase = false;
        }
      }
      else if(iLiability < TcutoffH) {
        if(nControl < ss) {
          iDH = 0;
          iIdx = ss + nControl;
          ++nControl;
        }
        else {
          rControl = false;
        }
      }
      
      // if there is case or control produced
      if (iIdx > -1) {
        gValue[iIdx] = ilg;
        eValue[iIdx] = ile;
        double imValue = scaleE * ile + R::rnorm(meanM, sdME);
        mValue[iIdx] = imValue;
        bValue[iIdx] = 1/(1+pow(2, -imValue));
        liability[iIdx] = iLiability;
        DH[iIdx] = iDH;
        // ++sCounter;
        // if(sCounter == sCounterTh) {
        //   sCounter = 0;
        //   Rcout << "\r"
        //         << std::fixed << std::setprecision(2)
        //         << "Getting sample (Case[" << (double (nCase) / ss * 100) << "%] : "
        //         << nCase << " of " << ss
        //         << ", Control[" << (double (nControl) / ss * 100) << "%]: "
        //         << nControl << " of " << ss << ")";
        // }
      }
      
    }
    // Rcout << "\r" << std::fixed << std::setprecision(2)
    //       << "Getting sample (Case[100.00%] : " << nCase << " of " << ss
    //       << ", Control[100.00%]: " << nControl << " of " << ss << ")";
    // Rcout << std::endl;
    
    return DataFrame::create(_["gValue"] = gValue,
                             _["eValue"] = eValue,
                             _["mValue"] = mValue,
                             _["bValue"] = bValue,
                             _["liability"] = liability,
                             _["DH"] = DH);
  }
  else {
    //** twin **//
    
    // twin methylation: multivariate normal distribution reference
    // and covariance matrix computation
    Function mvrnorm = Environment::namespace_env("MASS")["mvrnorm"];
    NumericVector muM(2, meanM);
    NumericMatrix covM(2, 2);
    covM(0, 0) = varME;
    covM(1, 1) = varME;
    covM(0, 1) = corME * varME;
    covM(1, 0) = corME * varME;
    
    
    // data.frame output
    NumericVector gValue(ss);
    NumericVector eValueD(ss);
    NumericVector mValueD(ss);
    NumericVector bValueD(ss);
    NumericVector liabilityD(ss);
    NumericVector DHD(ss);
    NumericVector eValueH(ss);
    NumericVector mValueH(ss);
    NumericVector bValueH(ss);
    NumericVector liabilityH(ss);
    NumericVector DHH(ss, 1.0);
    NumericVector deltaBvalue(ss);
    NumericVector deltaMvalue(ss);
    
    // number of discordant twins
    int nDT = 0;
    
    Rcout << "Getting twin sample" << std::endl;
    
    // // output counter
    // int sCounter = 0;
    // Rcout << "\r" << std::fixed << std::setprecision(2)
    //       << "Getting twin sample (discordant [0.00%] "
    //       << nDT << " of " << ss << ")";
    // // print when counter is
    // const int sCounterTh = 1000;
    
    while (nDT < ss) {
      
      double ilg = R::rnorm(0, sdG);
      
      double ile1 = R::rnorm(0, sdE);
      double ile2 = R::rnorm(0, sdE);
      
      double iLiability1 = ilg + ile1;
      double iLiability2 = ilg + ile2;
      
      if((iLiability1 > TcutoffD) & (iLiability2 < TcutoffH)) {
        gValue[nDT] = ilg;
        eValueD[nDT] = ile1;
        eValueH[nDT] = ile2;
        liabilityD[nDT] = iLiability1;
        liabilityH[nDT] = iLiability2;
        ++nDT;
        // ++sCounter;
        // if(sCounter == sCounterTh) {
        //   sCounter = 0;
        //   Rcout << "\r" << std::fixed << std::setprecision(2)
        //         << "Getting twin sample (discordant [" << (double (nDT) / ss * 100) << "%] "
        //         << nDT << " of " << ss << ")";
        // }
      }
      else if((iLiability1 < TcutoffH) & (iLiability2 > TcutoffD)) {
        gValue[nDT] = ilg;
        eValueD[nDT] = ile2;
        eValueH[nDT] = ile1;
        liabilityD[nDT] = iLiability2;
        liabilityH[nDT] = iLiability1;
        ++nDT;
        // ++sCounter;
        // if(sCounter == sCounterTh) {
        //   sCounter = 0;
        //   Rcout << "\r" << std::fixed << std::setprecision(2)
        //         << "Getting twin sample (discordant [" << (double (nDT) / ss * 100) << "%] "
        //         << nDT << " of " << ss << ")";
        // }
      }
    }
    
    NumericMatrix mErr = mvrnorm(ss, muM, covM);
    
    for(int i = 0; i < ss; ++i) {
      double imValueD = scaleE * eValueD[i] + mErr(i, 0);
      double imValueH = scaleE * eValueH[i] + mErr(i, 1);
      double ibValueD = 1/(1+pow(2, -imValueD));
      double ibValueH = 1/(1+pow(2, -imValueH));
      mValueD[i] = imValueD;
      mValueH[i] = imValueH;
      bValueD[i] = ibValueD;
      bValueH[i] = ibValueH;
      deltaMvalue[i] = imValueD - imValueH;
      deltaBvalue[i] = ibValueD - ibValueH;
    }
    
    // Rcout << "\r" << std::fixed << std::setprecision(2)
    //       << "Getting twin sample (discordant [100.00%] "
    //       << nDT << " of " << ss << ")";
    // Rcout << std::endl;
    
    return DataFrame::create(
      _["gValue"] = gValue,
      _["eValueD"] = eValueD,
      _["eValueH"] = eValueH,
      _["mValueD"] = mValueD,
      _["mValueH"] = mValueH,
      _["bValueD"] = bValueD,
      _["bValueH"] = bValueH,
      _["liabilityD"] = liabilityD,
      _["liabilityH"] = liabilityH,
      _["DHD"] = DHD,
      _["DHH"] = DHH,
      _["deltaBvalue"] = deltaBvalue,
      _["deltaMvalue"] = deltaMvalue);
  }
}
