#include <Rcpp.h>
using namespace Rcpp;

// compute median ignoring missing values
double median(NumericVector x) {
	NumericVector y = na_omit(x);
	int pos = y.size()/2;
	std::nth_element(y.begin(), y.begin()+pos, y.end());
	return y[pos];
}

// ancient compilers (eg. on OS X < 10.7) may not have std::log2(), so...
double log2(double x) {
	return log(x)/log(2.0);
}

// [[Rcpp::export]]
List tQN_Cpp( NumericVector r, NumericVector theta, 
NumericVector A_theta, NumericVector A_r,
NumericVector B_theta, NumericVector B_r, 
NumericVector H_theta, NumericVector H_r) {

	int n = r.size();
	NumericVector baf = NumericVector(n);
	NumericVector lrr = NumericVector(n);
	
	// get median thetas for AA,BB clusters
	double med_tAA = median(A_theta);
	double med_tBB = median(B_theta);
	
	// check existence of reference thetas
	LogicalVector e_tAA = !is_na(A_theta);
	LogicalVector e_tAB = !is_na(H_theta);
	LogicalVector e_tBB = !is_na(B_theta);
	
	// check existence of sample thetas
	LogicalVector theta_na = is_na(theta);
	
	// loop on markers
	for (int i = 0; i < n; i++) {
		
		// theta doesn't exist? set everything missing and move on
		if (theta_na[i]) {
			baf[i] = NAN;
			lrr[i] = NAN;
			continue;
		}
		
		// make convenience copies for this loop iteration; should be cheap
		double th = theta[i];
		double rr = r[i];
		double rAA = A_r[i];
		double rAB = H_r[i];
		double rBB = B_r[i];
		double tAA = A_theta[i];
		double tAB = H_theta[i];
		double tBB = B_theta[i];
		
		// check for A/B allele swaps
		if (e_tAA[i] && e_tBB[i] && tAA > tBB) {
			double r_tmp = rBB;
			double t_tmp = tBB;
			rBB = rAA;
			tBB = tAA;
			rAA = r_tmp;
			tAA = t_tmp;
		}
		
		// 0: Test for inconsistencies between tAA/tAB/tBB and rAA/rAB/rBB
		if (((e_tAA[i] && e_tAB[i]) && tAA > tAB) ||
				((e_tAA[i] && e_tBB[i]) && tAA > tBB) ||
				((e_tAB[i] && e_tBB[i]) && tAB > tBB)) {
			baf[i] = NAN;
			lrr[i] = NAN;
		}
		// 1: Triple blank SNP
		else if (!(e_tAA[i] || e_tAB[i] || e_tBB[i])) {
			baf[i] = NAN;
			lrr[i] = NAN;
		}
		// 2: Blank for AB, AA, while positive for BB
		else if (!(e_tAA[i] || e_tAB[i]) && e_tBB[i]) {
			if (th >= tBB) {
				baf[i] = 1.0;
				lrr[i] = rBB <= 0 ? NAN : log2(rr / rBB);
			}
			else {
				baf[i] = NAN;
				lrr[i] = NAN;
			}
		}
		// 3: Blank for AB, BB, while positive for AA
		else if (e_tAA[i] && !(e_tAB[i] || e_tBB[i])) {
			if (th <= tAA) {
				baf[i] = 0.0;
				lrr[i] = tAA <= 0 ? NAN : log2(rr / tAA);
			}
			else {
				baf[i] = NAN;
				lrr[i] = NAN;
			}
		}
		// 4: Blank for AB while positive for AA & BB
		else if (e_tAA[i] && !e_tAB[i] && e_tBB[i]) {
			// No AB cluster exist for this SNP, while AA & BB exists.
			// Set it to the closest of AA or BB
			int min_index = 0;
			if (std::abs(tBB - th) < std::abs(tAA - th)) {
				min_index = 1;
			}
			if (min_index == 1 && th < tAA) {
				baf[i] = 0.0;
				lrr[i] = tAA <= 0 ? NAN : log2(rr / tAA);
			}
			else if (min_index != 1 && th >= tBB) {
				baf[i] = 1.0;
				lrr[i] = rBB <= 0 ? NAN : log2(rr / rBB);
			}
			else {
				baf[i] = NAN;
				lrr[i] = NAN;
			}
		}
		// 5: Blank for AA while positive for AB & BB
		else if (!e_tAA[i] && e_tAB[i] && e_tBB[i]) {
			if (th >= tBB) {
				baf[i] = 1.0;
				lrr[i] = NAN;
			}
			// 5.1: SNP is "correctly between" ref$AB_T_Mean and ref$BB_T_Mean
			else if (th >= tAB) {
				// interpolate as SNP is expected to be between ref$AB_T_Mean and ref$BB_T_Mean
				baf[i] = 0.5 + 0.5 * (th - tAB) / (tBB - tAB);
				double eR = rAB + ((th - tAB) * (rBB - rAB) / (tBB - tAB));
				lrr[i] = eR <= 0 ? NAN : log2(rr / eR);
			}
			// 5.2: Heterozygous SNP is subjected to deletion or UPD of allele B making it unexectedly to be 
			// between ref$AA_T_Mean and ref$AB_T_Mean where it normally should not NOT BE.
			else {
				baf[i] = th < med_tAA ? 0.0 : 0.5 * (th - med_tAA) / (tAB - med_tAA);
				lrr[i] = NAN;
			}
		}
		// 6: Blank for BB while positive for AA & AB
		else if (e_tAA[i] && e_tAB[i] && !e_tBB[i]) {
			if (th < tAA) {
				baf[i] = 0.0;
				lrr[i] = NAN;
			}
			// 6.1: SNP is "correctly between" ref$AA_T_Mean and ref$AB_T_Mean
			else if (th <= tAB) {
				// interpolate as SNP is expected to be between ref$AB_T_Mean and ref$BB_T_Mean
				baf[i] = 0.5* (th - tAA) / (tAB - tAA);
				double eR = rAA + ((th - tAA) * (rAB - rAA) / (tAB - tAA));
				lrr[i] = eR <= 0 ? NAN : log2(rr / eR);
			}
			// 2: Heterozygous SNP is subjected to deletion or UPD of allele A making it unexectedly to be 
			// between ref$AB_T_Mean and ref$BB_T_Mean where it normally should not NOT BE.
			else {
				baf[i] = th > med_tBB ? 1.0 : 0.5 + 0.5 * (th - tAB) / (med_tBB - tAB);
				lrr[i] = NAN;
			}
		}
		// 7: positive for AA & BB & AB, Illumina style calculation
		else {
			double eR = 0.0;
			if (th < tAB) {
				baf[i] = th < tAA ? 0.0 : 0.5 * (th - tAA) / (tAB - tAA);
				eR = rAA + ((th - tAA) * (rAB - rAA) / (tAB - tAA));
			}
			else {
				baf[i] = th >= tBB ? 1.0 : 0.5 + 0.5 * (th - tAB) / (tBB - tAB);
				eR = rAB + ((th - tAB) * (rBB - rAB) / (tBB - tAB));
			}
			lrr[i] = eR <= 0 ? NAN : log2(rr / eR);
		}
	}
	
	//if (adj_lrr) {
	//	lrr = lrr / R_mean;
	//}
	
	// wrap it in a list and return
	List rez;
	rez["LRR"] = lrr;
	rez["BAF"] = baf;
	return rez;

} // end tQN_Cpp()
