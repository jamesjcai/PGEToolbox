/*
 * vcf_file_output.cpp
 *
 *  Created on: Aug 28, 2009
 *      Author: Adam Auton
 *      ($Revision: 119 $)
 */
#include "vcf_file.h"


// Output as PLINK formatted PED/MAP files.
void vcf_file::output_as_plink(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output as PLINK.");

	printLOG("Writing PLINK file ... ");
	string ped_file = output_file_prefix + ".ped";
	string map_file = output_file_prefix + ".map";

	ofstream PED(ped_file.c_str());
	if (!PED.is_open()) error("Could not open output file: " + ped_file, 12);

	unsigned int ui, uk, s;
	vector<string> alleles;
	pair<int, int> genotype;
	string vcf_line;
	vcf_entry e(N_indv);
	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		PED << indv[ui] << "\t" << indv[ui] << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0;
		uk = 2*ui;
		for (s=0; s<N_entries; s++)
		{
			if (include_entry[s] == false)
				continue;

			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_basic_entry(true);

			if (e.get_N_alleles() == 2)	// Only output sites with one alternative allele
			{
				e.get_alleles_vector(alleles);
				genotype = make_pair(-1,-1);
				if (include_genotype[s][ui] == true)
				{
					e.parse_genotype_entry(ui, true);
					e.get_indv_GENOTYPE_ids(ui, genotype);
				}

				if (genotype.first == -1)
					PED << "\t0";
				else
					PED << "\t" << alleles[genotype.first];

				if (genotype.second == -1)
				{
					if (e.get_indv_PHASE(ui) == '/')
						PED << "\t0";
					else
						PED << "\t" << alleles[genotype.first];	// Male X-chr, Y-chr etc
				}
				else
					PED << "\t" << alleles[genotype.second];
			}
		}
		PED << endl;
	}

	PED.close();

	ofstream MAP(map_file.c_str());
	if (!MAP.is_open()) error("Could not open output file: " + map_file, 12);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		if (e.get_N_alleles() == 2)	// Only output sites with one alternative allele
		{
			if (e.get_ID() == ".")
				MAP << e.get_CHROM() << "\t" << e.get_POS() << "\t0\t" << e.get_POS() << endl;
			else
				MAP << e.get_CHROM() << "\t" << e.get_ID() << "\t0\t" << e.get_POS() << endl;
		}
	}

	MAP.close();
	printLOG("Done.\n");
}

// Output as a simple 0/1/2 matrix
void vcf_file::output_as_012_matrix(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output as 0/1/2 matrix.");

	printLOG("Writing 012 matrix file ... ");
	string ped_file = output_file_prefix + ".012";
	string map_file = output_file_prefix + ".012.pos";
	string fam_file = output_file_prefix + ".012.indv";

	ofstream PED(ped_file.c_str());
	if (!PED.is_open()) error("Could not open output file: " + ped_file, 12);
	string allele1, allele2;

	ofstream FAM(fam_file.c_str());
	if (!FAM.is_open()) error("Could not open output file: " + fam_file, 12);

	unsigned int ui, s, uk;
	pair<int, int> genotype;
	string vcf_line;
	vcf_entry e(N_indv);
	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		FAM << indv[ui] << endl;
		PED << ui;
		uk = 2*ui;
		for (s=0; s<N_entries; s++)
		{
			if (include_entry[s] == false)
				continue;

			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_basic_entry(true);

			if (e.get_N_alleles() == 2)	// Only output sites with one alternative allele
			{
				genotype = make_pair(-1,-1);
				if (include_genotype[s][ui] == true)
				{
					e.parse_genotype_entry(ui, true);
					e.get_indv_GENOTYPE_ids(ui, genotype);
				}

				if ((genotype.first == -1) && (genotype.second == -1))
					PED << "\t-1";	// Missing data
				else if ((genotype.first == 0) && (genotype.second == 0))
					PED << "\t0";	// No copies of the alternative allele
				else
				{
					if ((genotype.first == 1) && (genotype.second == 1))
						PED << "\t2";	// Two copies of the alternative allele
					else
						PED << "\t1";	// Must be one copy of the alternative allele.
				}
			}
		}
		PED << endl;
	}

	FAM.close();
	PED.close();

	ofstream MAP(map_file.c_str());
	if (!MAP.is_open()) error("Could not open output file: " + map_file, 12);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		if (e.get_N_alleles() == 2)	// Only output sites with one alternative allele
		{
			MAP << e.get_CHROM() << "\t" << e.get_POS() << endl;
		}
	}

	MAP.close();
	printLOG("Done.\n");
}

// Output statistics of frequency at each site
void vcf_file::output_frequency(const string &output_file_prefix, bool output_counts)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output Frequency Statistics.");

	printLOG("Outputting Frequency Statistics...\n");
	string output_file = output_file_prefix + ".frq";
	if (output_counts)
		output_file += ".count";

	ofstream out(output_file.c_str());
	if (!out.is_open()) error("Could not open output file: " + output_file, 12);
	out << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:";
	if (output_counts)
		out << "COUNT}" << endl;
	else
		out << "FREQ}" << endl;

	unsigned int ui;
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned int N_alleles;
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		e.parse_genotype_entries(true);
		N_alleles = e.get_N_alleles();

		e.get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);

		out << e.get_CHROM() << "\t" << e.get_POS() << "\t" << N_alleles << "\t" << N_non_missing_chr;
		if (output_counts)
		{
			out << "\t" << e.get_REF() << ":" << allele_counts[0];
			for (ui=1; ui<N_alleles; ui++)
			{
				out << "\t" << e.get_ALT_allele(ui-1) << ":" << allele_counts[ui];
			}
			out << endl;
		}
		else
		{
			double freq;
			freq = allele_counts[0] / (double)N_non_missing_chr;
			out << "\t" << e.get_REF() << ":" << freq;
			for (ui=1; ui<N_alleles; ui++)
			{
				freq = allele_counts[ui] / (double)N_non_missing_chr;
				out << "\t" << e.get_ALT_allele(ui-1) << ":" << freq;
			}
			out << endl;
		}
	}
	out.close();
}

// Output statistics on Heterozygosity for each individual
void vcf_file::output_het(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output Heterozygosity Statistics.");
	// Following the calculations in PLINK....
	// Note this assumes Biallelic SNPs.

	printLOG("Outputting Individual Heterozygosity (but only for biallelic loci)\n");

	string output_file = output_file_prefix + ".het";
	ofstream out(output_file.c_str());
	if (!out.is_open()) error("Could not open output file: " + output_file, 12);
	out << "INDV\tO(HOM)\tE(HOM)\tN(NM)\tF" << endl;

	// P(Homo) = F + (1-F)P(Homo by chance)
	// P(Homo by chance) = p^2+q^2 for a biallelic locus.
	// For an individual with N genotyped loci, we
	//   1. count the total observed number of loci which are homozygous (O),
	//   2. calculate the total expected number of loci homozygous by chance (E)
	// Then, using the method of moments, we have
	//    O = NF + (1-F)E
	// Which rearranges to give
	//    F = (O-E)/(N-E)

	// First, calc frequency of each site (should really move this to a subroutine)
	vector<double> freq(N_entries, 0.0);
	unsigned int ui, s;
	vector<int> allele_counts;
	vector<unsigned int> N_non_missing_chr(N_entries,0);
	string vcf_line;
	vcf_entry e(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
			continue;

		e.parse_genotype_entries(true);
		// Frequency of non-reference allele
		e.get_allele_counts(allele_counts, N_non_missing_chr[s], include_indv, include_genotype[s]);

		if (N_non_missing_chr[s] > 0)
			freq[s] = allele_counts[1] / double(N_non_missing_chr[s]);
		else
			freq[s] = -1;
	}

	vector<int> N_sites_included(N_indv, 0);
	vector<int> N_obs_hom(N_indv, 0);
	vector<double> N_expected_hom(N_indv, 0.0);
	pair<int, int> alleles;

	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
			continue;

		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if ((freq[s] <= numeric_limits<double>::epsilon())  || (1.0 - freq[s] <= numeric_limits<double>::epsilon()))
				continue;

			e.parse_genotype_entry(ui, true);
			if (include_genotype[s][ui] == true)
			{
				e.get_indv_GENOTYPE_ids(ui, alleles);
				if ((alleles.first != -1) && (alleles.second != -1))
				{
					N_sites_included[ui]++;
					if (alleles.first == alleles.second)
						N_obs_hom[ui]++;
				}

				/////////////////////////
				// Expected homozygosity
				// E = 1 - (2pq . 2N/(2N-1))
				// (Using Nei's unbiased estimator)
				N_expected_hom[ui] += 1.0 - (2.0 * freq[s] * (1.0 - freq[s]) * (N_non_missing_chr[s] / (N_non_missing_chr[s] - 1.0)));
			}
		}
	}

	out.setf(ios::fixed,ios::floatfield);
	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		if (N_sites_included[ui] > 0)
		{
			double F = (N_obs_hom[ui] - N_expected_hom[ui]) / double(N_sites_included[ui] - N_expected_hom[ui]);
			out << indv[ui] << "\t" << N_obs_hom[ui] << "\t";
			out.precision(1);
			out << N_expected_hom[ui] << "\t";
			out.precision(5);
			out << N_sites_included[ui] << "\t" << F << endl;
		}
	}

	out.close();
}

// Output HWE statistics for each site
void vcf_file::output_hwe(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output HWE Statistics.");
	// Note this assumes Biallelic SNPs.
	printLOG("Outputting HWE statistics (but only for biallelic loci)\n");

	string output_file = output_file_prefix + ".hwe";
	ofstream out(output_file.c_str());
	if (!out.is_open()) error("Could not open output file: " + output_file, 12);
	out << "CHR\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tChiSq\tP" << endl;

	/* PLINK code:
	// b11 = Nhom1, b12 = Nhet, b22 = Nhom2
	double tot = b11 + b12 + b22;
	double exp_11 = freq * freq * tot;
	double exp_12 = 2 * freq * (1-freq) * tot;
	double exp_22 = (1-freq) * (1-freq) * tot;

	double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
			    + ( (b12-exp_12)*(b12-exp_12) ) / exp_12
			    + ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;

	p = chiprobP(chisq,1);
	*/

	double freq;
	unsigned int b11, b12, b22;
	double exp_11, exp_12, exp_22;
	unsigned int s;
	double chisq;
	double tot;
	double p;
	unsigned int precision = out.precision();
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	string vcf_line;
	vcf_entry e(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() > 2)
			continue;

		e.parse_genotype_entries(true);

		e.get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
		freq = allele_counts[0] / (double)N_non_missing_chr;
		e.get_genotype_counts(include_indv, include_genotype[s], b11, b12, b22);
		tot = b11 + b12 + b22;
		exp_11 = freq * freq * tot;
		exp_12 = 2.0 * freq * (1.0-freq) * tot;
		exp_22 = (1.0-freq) * (1.0-freq) * tot;

		chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
				+ ( (b12-exp_12)*(b12-exp_12) ) / exp_12
				+ ( (b22-exp_22)*(b22-exp_22) ) / exp_22;

		p = vcf_entry::SNPHWE(b12, b11, b22);
		out << e.get_CHROM() << "\t" << e.get_POS();
		out << "\t" << b11 << "/" << b12 << "/" << b22;
		out.precision(2);
		out << fixed << "\t" << exp_11 << "/" << exp_12 << "/" << exp_22;
		out.precision(precision);
		out << "\t" << chisq << "\t" << p << endl;
	}
}

// Output information regarding the mean depth for each individual
void vcf_file::output_individuals_by_mean_depth(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output Individuals by Mean Depth Statistics.");

	printLOG("Outputting Mean Depth by Individual\n");
	string output = output_file_prefix + ".idepth";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open Individual Depth Output File: " + output, 2);
	out << "INDV\tN_SITES\tMEAN_DEPTH" << endl;
	unsigned int ui;
	vector<double> depth_sum(N_indv, 0.0);
	vector<int> count(N_indv, 0);
	int depth;
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);

		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == true)
			{
				e.parse_genotype_entry(ui, false, false, true);
				depth = e.get_indv_DEPTH(ui);
				if (depth >= 0)
				{
					depth_sum[ui] += depth;
					count[ui]++;
				}
			}
		}
	}

	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		double mean_depth = depth_sum[ui] / count[ui];
		out << indv[ui] << "\t" << count[ui] << "\t" << mean_depth << endl;
	}

	out.close();
}

// Output as IMPUTE format
void vcf_file::output_as_IMPUTE(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output IMPUTE format.");

	printLOG("Outputting in IMPUTE format (bi-allelic, completely phased SNPs only)\n");
	unsigned int s, ui;
	string legend_file = output_file_prefix + ".impute.legend";
	string haplotype_file = output_file_prefix + ".impute.hap";
	string indv_file = output_file_prefix + ".impute.hap.indv";
	ofstream legend(legend_file.c_str());
	if (!legend.is_open())
		error("Could not open IMPUTE Legend Output File: " + legend_file, 2);
	legend << "ID pos allele0 allele1" << endl;

	ofstream hap(haplotype_file.c_str());
	if (!hap.is_open())
		error("Could not open IMPUTE Haplotype Output File: " + haplotype_file, 2);

	ofstream indv_out(indv_file.c_str());
	if (!indv_out.is_open())
		error("Could not open IMPUTE Individual Output File: " + indv_file, 2);

	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		indv_out << indv[ui] << endl;
	}
	indv_out.close();

	pair<int, int> alleles;
	string vcf_line;
	vcf_entry e(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
			continue;

		// Exclude entries with missing data and/or unphased
		bool missing = false;
		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == false)
			{
				missing = true;
				break;
			}

			e.parse_genotype_entry(ui, true);
			e.get_indv_GENOTYPE_ids(ui, alleles);
			if ((alleles.first == -1) || (alleles.second == -1))
			{
				missing = true;
				break;
			}

			if (e.get_indv_PHASE(ui) != '|')
			{
				missing = true;
				break;
			}
		}
		if (missing == true)
			continue;

		if (e.get_ID() == ".")
		{
			legend << e.get_CHROM() << "-" << e.get_POS() << " " << e.get_POS() << " " << e.get_REF() << " " << e.get_ALT_allele(0) << endl;
		}
		else
			legend << e.get_ID() << " " << e.get_POS() << " " << e.get_REF() << " " << e.get_ALT_allele(0) << endl;

		bool first = true;
		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.parse_genotype_entry(ui, true);
			e.get_indv_GENOTYPE_ids(ui, alleles);
			if (first == true)
			{
				hap << alleles.first << " " << alleles.second;
				first = false;
			}
			else
				hap << " " << alleles.first << " " << alleles.second;
		}
		hap << endl;
	}

	hap.close();
	legend.close();
}

// Output as LDhat format
void vcf_file::output_as_LDhat(const string &output_file_prefix, const string &chr)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output LDhat format.");

	printLOG("Outputting in LDhat format\n");
	if (chr == "")
		error("Require chromosome for LDhat output", 10);

	string locs_file = output_file_prefix + ".ldhat.locs";
	string sites_file = output_file_prefix + ".ldhat.sites";

	ofstream locs(locs_file.c_str());
	if (!locs.is_open())
		error("Could not open LDhat locs Output File: " + locs_file, 2);

	ofstream sites(sites_file.c_str());
	if (!sites.is_open())
		error("Could not open LDhat sites Output File: " + sites_file, 2);

	unsigned int s, ui, k;
	int max_pos = -1;
	unsigned int n_sites=0;
	unsigned int n_indv = N_kept_individuals();

	vcf_entry e(N_indv);
	string vcf_line;
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry();

		max_pos = max(e.get_POS(), max_pos);
		n_sites++;
	}

	locs << n_sites;
	locs.setf(ios::fixed,ios::floatfield);
	locs.precision(4);
	locs << "\t" << max_pos / 1000.0 << "\tL" << endl;
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry();

		locs << e.get_POS() / 1000.0 << endl;
	}

	sites << n_indv << "\t" << n_sites << "\t1" << endl;
	// TODO: This is quite painfully slow as the get_vcf_entry is in the inner loop. Can I speed it up?
	pair<int, int> alleles;
	unsigned int max_k = 2;
	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		// Check to see if this is a male X-chr.
		for (s=0; s<N_entries; s++)
		{
			if (include_entry[s] == false)
				continue;

			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_genotype_entry(ui, true);

			e.get_indv_GENOTYPE_ids(ui, alleles);

			if ((alleles.first != -1) && (alleles.second != -1))
			{	// Diploid
				max_k = 2;
				break;
			}

			if ((alleles.first != -1) && (alleles.second == -1) && (e.get_indv_PHASE(ui) == '|'))
			{	// Haploid
				max_k = 1;
				break;
			}
		}

		for (k=0; k<max_k; k++)
		{
			sites << ">" << indv[ui] << "-" << k+1 << endl;
			for (s=0; s<N_entries; s++)
			{
				if (include_entry[s] == false)
					continue;

				get_vcf_entry(s, vcf_line);
				e.reset(vcf_line);

				e.parse_genotype_entry(ui, true);

				e.get_indv_GENOTYPE_ids(ui, alleles);
				int geno;
				if (k == 0)
					geno = alleles.first;
				else
					geno = alleles.second;

				if (geno != -1)
					sites << geno;
				else
					sites << "?";
			}
			sites << endl;
		}
	}

	sites.close();
	locs.close();
}

// Output SNP density
void vcf_file::output_SNP_density(const string &output_file_prefix, int bin_size)
{
	if (bin_size <= 0)
		return;
	printLOG("Outputting SNP density\n");

	string output = output_file_prefix + ".snpden";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open SNP Density Output File: " + output, 2);

	// Find maximum position
	unsigned int s;
	map<string, int> max_pos;
	string vcf_line;
	string CHROM;
	vcf_entry e(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_basic_entry();

			CHROM = e.get_CHROM();

			if (max_pos.find(CHROM) != max_pos.end())
			{
				if (e.get_POS() > max_pos[CHROM])
					max_pos[CHROM] = e.get_POS();
			}
			else
				max_pos[CHROM] = e.get_POS();
		}
	}

	map<string, int>::iterator it;

	unsigned int N_bins;
	map<string, vector<int> > bins;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		N_bins = (unsigned int)((max_pos[CHROM] + bin_size) / double(bin_size));
		bins[CHROM].resize(N_bins, 0);
	}


	unsigned int idx;
	double C = 1.0 / double(bin_size);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_basic_entry();

			CHROM = e.get_CHROM();

			idx = (unsigned int)(e.get_POS() * C);
			bins[CHROM][idx]++;
		}
	}

	out << "CHROM\tBinStart\tSNP_count\tSNPs/kb" << endl;
	double sum1=0.0, sum2=0.0;
	int bin_tot;
	C = 1000.0 / bin_size;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		sum2 += max_pos[CHROM];
		for (s=0; s<bins[CHROM].size(); s++)
		{
			bin_tot = bins[CHROM][s];
			sum1 += bin_tot;
			out << CHROM << "\t" << s*bin_size << "\t" << bin_tot << "\t" << bin_tot * C << endl;
		}
	}
	out.close();

	double mean_SNP_density = sum1 / sum2 * 1000;
	printLOG("Mean SNP density: " + dbl2str(mean_SNP_density, 5) + " SNPs / kb\n");

}

// Output missingness by individual and site
void vcf_file::output_missingness(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output Missingness Statistics.");

	printLOG("Outputting Individual Missingness\n");
	string output = output_file_prefix + ".imiss";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open Individual Missingness Output File: " + output, 3);

	out << "INDV\tN_DATA\tN_SITES_FILTERED\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS" << endl;
	unsigned int ui, s;
	unsigned int N_missing, N_tot;
	unsigned int N_site_filtered, N_geno_filtered;
	pair<int, int> alleles;
	string vcf_line;
	vcf_entry e(N_indv);
	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		N_missing = 0;
		N_tot = 0;
		N_site_filtered = 0;
		N_geno_filtered = 0;
		for (s=0; s<N_entries; s++)
		{
			if (include_entry[s] == false)
			{
				N_site_filtered++;
				continue;
			}
			if (include_genotype[s][ui] == false)
			{
				N_geno_filtered++;
				continue;
			}

			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_genotype_entry(ui, true);
			e.get_indv_GENOTYPE_ids(ui, alleles);
			if (alleles.first == -1)
				N_missing++;
			N_tot+=1;
		}
		out << indv[ui] << "\t" << N_tot << "\t" << N_site_filtered << "\t";
		out << N_geno_filtered << "\t" << N_missing << "\t" << N_missing / double(N_tot) << endl;
	}
	out.close();

	printLOG("Outputting Site Missingness\n");
	output = output_file_prefix + ".lmiss";
	out.open(output.c_str());
	if (!out.is_open())
		error("Could not open Site Missingness Output File: " + output, 4);

	out << "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS" << endl;
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry();

		N_missing = 0;
		N_tot = 0;
		N_geno_filtered = 0;
		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			if (include_genotype[s][ui] == false)
			{
				N_geno_filtered++;
				continue;
			}

			e.parse_genotype_entry(ui, true);
			e.get_indv_GENOTYPE_ids(ui, alleles);
			if (alleles.first == -1)
			{
				N_missing++;
			}

			if (alleles.second == -1)
			{
				N_missing++;
			}
			N_tot+=2;

			if ((alleles.second == -1) && (e.get_indv_PHASE(ui) == '|'))
			{	// Phased missing genotypes indicate haploid genome
				N_tot--;
			}
		}
		out << e.get_CHROM() << "\t" << e.get_POS() << "\t" << N_tot << "\t" << N_geno_filtered << "\t";
		out << N_missing << "\t" << double(N_missing) / double(N_tot) << endl;
	}

	out.close();
}

// Output pairwise LD statistics, using traditional r^2. Requires phased haplotypes.
void vcf_file::output_haplotype_r2(const string &output_file_prefix, int snp_window_size, int bp_window_size, double min_r2)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output LD Statistics.");

	unsigned int s, s2;
	unsigned int ui;

	printLOG("Outputting Pairwise LD (phased bi-allelic only)\n");
	string output = output_file_prefix + ".hap.ld";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open LD Output File: " + output, 3);

	out << "CHR\tPOS1\tPOS2\tN_CHR\tR^2" << endl;

	unsigned int chr_count;
	double r2;
	int sx, sy;
	double X, X2, Y, Y2, XY;
	double var1, var2, cov12;
	pair<int,int> geno1, geno2;
	string vcf_line, vcf_line2;
	vcf_entry e(N_indv), e2(N_indv);
	for (s=0; s<(N_entries-1); s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
			continue;

		for (s2 = s+1; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;

			if (int(s2 - s) > snp_window_size)
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			get_vcf_entry(s2, vcf_line2);
			e2.reset(vcf_line2);
			e2.parse_basic_entry(true);

			if (e.get_CHROM() != e2.get_CHROM())
			{
				s2 = N_entries;	// No need to go any further (assuming SNPs are sorted)
				continue;
			}

			if ((e2.get_POS() - e.get_POS()) > bp_window_size)
			{
				s2 = N_entries;	// No need to go any further (assuming SNPs are sorted)
				continue;
			}

			if (e2.get_N_alleles() != 2)
				continue;

			X=0, X2=0; Y=0; Y2=0; XY=0;
			chr_count = 0;
			for (ui=0; ui<N_indv; ui++)
			{
				if ((include_indv[ui] == false) || (include_genotype[s][ui] == false) || (include_genotype[s2][ui] == false))
					continue;

				e.parse_genotype_entry(ui, true);
				e.get_indv_GENOTYPE_ids(ui, geno1);

				e2.parse_genotype_entry(ui, true);
				e2.get_indv_GENOTYPE_ids(ui, geno2);

				if ((e.get_indv_PHASE(ui) != '|') || (e2.get_indv_PHASE(ui) != '|'))
					error("Require phased haplotypes for r^2 calculation (use --phased)\n");

				for (unsigned int c=0; c<2; c++)
				{
					int allele1, allele2;
					if (c==0)
					{
						allele1 = geno1.first;
						allele2 = geno2.first;
					}
					else
					{
						allele1 = geno1.second;
						allele2 = geno2.second;
					}

					if ((allele1 == -1) || (allele2 == -1))
						continue;

					sx=0, sy=0;
					if (allele1 == 0)
						sx += 1;

					if (allele2 == 0)
						sy += 1;

					X += sx;
					Y += sy;
					XY += sx*sy;

					sx *= sx;
					sy *= sy;

					X2 += sx;
					Y2 += sy;

					chr_count++;
				}
			}

			X /= chr_count;
			X2 /= chr_count;
			Y /= chr_count;
			Y2 /= chr_count;
			XY /= chr_count;

			var1 = X2 - X*X;
			var2 = Y2 - Y*Y;
			cov12 = XY - X*Y;

			r2 = cov12 * cov12 / (var1 * var2);

			if (r2 < min_r2)
				continue;

			out << e.get_CHROM() << "\t" << e.get_POS() << "\t" << e2.get_POS() << "\t" << chr_count << "\t" << r2 << endl;
		}
	}
	out.close();
}

// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
// correlation coefficient between genotypes numbered as 0, 1, 2.
void vcf_file::output_genotype_r2(const string &output_file_prefix, int snp_window_size, int bp_window_size, double min_r2)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output LD Statistics.");

	unsigned int s, s2;
	unsigned int ui;

	printLOG("Outputting Pairwise LD (bi-allelic only)\n");
	string output = output_file_prefix + ".geno.ld";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open LD Output File: " + output, 3);

	out << "CHR\tPOS1\tPOS2\tN_INDV\tR^2" << endl;

	unsigned int indv_count;
	double r2;
	int sx, sy;
	double X, X2, Y, Y2, XY;
	double var1, var2, cov12;
	pair<int,int> geno1, geno2;
	string vcf_line, vcf_line2;
	vcf_entry e(N_indv), e2(N_indv);
	for (s=0; s<(N_entries-1); s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry();

		for (s2 = s+1; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;

			if (int(s2 - s) > snp_window_size)
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			get_vcf_entry(s2, vcf_line2);
			e2.reset(vcf_line2);
			e2.parse_basic_entry();

			if (e.get_CHROM() != e2.get_CHROM())
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			if ((e2.get_POS() - e.get_POS()) > bp_window_size)
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			X=0, X2=0; Y=0; Y2=0; XY=0;
			indv_count = 0;
			for (ui=0; ui<N_indv; ui++)
			{
				if ((include_indv[ui] == false) || (include_genotype[s][ui] == false) || (include_genotype[s2][ui] == false))
					continue;

				e.parse_genotype_entry(ui, true);
				e.get_indv_GENOTYPE_ids(ui, geno1);
				if ((geno1.first == -1) || (geno1.second == -1))
					continue;

				e2.parse_genotype_entry(ui, true);
				e2.get_indv_GENOTYPE_ids(ui, geno2);
				if ((geno2.first == -1) || (geno2.second == -1))
					continue;

				sx=0, sy=0;
				if (geno1.first == geno1.second)
				{
					if (geno1.first == 0)
					{
						sx = 2;
					}
				}
				else
					sx = 1;

				if (geno2.first == geno2.second)
				{
					if (geno2.first == 0)
					{
						sy = 2;
					}
				}
				else
					sy = 1;

				X += sx;
				Y += sy;
				XY += sx*sy;

				sx *= sx;
				sy *= sy;

				X2 += sx;
				Y2 += sy;

				indv_count++;
			}

			X /= indv_count;
			X2 /= indv_count;
			Y /= indv_count;
			Y2 /= indv_count;
			XY /= indv_count;

			var1 = X2 - X*X;
			var2 = Y2 - Y*Y;
			cov12 = XY - X*Y;

			r2 = cov12 * cov12 / (var1 * var2);

			if (r2 < min_r2)
				continue;

			out << e.get_CHROM() << "\t" << e.get_POS() << "\t" << e2.get_POS() << "\t" << indv_count << "\t" << r2 << endl;
		}
	}
	out.close();
}

void vcf_file::output_singletons(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output Singletons.");

	printLOG("Outputting Singleton Locations\n");
	string output = output_file_prefix + ".singletons";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open Singleton Output File: " + output, 3);

	out << "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV" << endl;

	unsigned int ui;
	int a;
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned int N_alleles;
	pair<int, int> geno;
	string allele;
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		e.parse_genotype_entries(true);

		e.get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
		N_alleles = e.get_N_alleles();

		for (a=0; a<(signed)N_alleles; a++)
		{
			if (allele_counts[a] == 1)
			{	// Singleton
				for (ui=0; ui<N_indv; ui++)
				{
					if (include_indv[ui] == false)
						continue;
					e.get_indv_GENOTYPE_ids(ui, geno);
					if ((geno.first == a) || (geno.second == a))
					{
						e.get_allele(a, allele);
						out << e.get_CHROM() << "\t" << e.get_POS() << "\tS\t" << allele << "\t" << indv[ui] << endl;
						ui=N_indv;
						break;
					}
				}
			}
			else if (allele_counts[a] == 2)
			{	// Possible doubleton
				for (ui=0; ui<N_indv; ui++)
				{
					if (include_indv[ui] == false)
						continue;
					e.get_indv_GENOTYPE_ids(ui, geno);
					if ((geno.first == a) || (geno.second == a))
					{
						e.get_allele(a, allele);
						out << e.get_CHROM() << "\t" << e.get_POS() << "\tD\t" << allele << "\t" << indv[ui] << endl;
						ui=N_indv;
						break;
					}
				}
			}
		}
	}

	out.close();
}

void vcf_file::output_INFO_for_each_site(const string &output_file_prefix, const vector<string> &INFO_to_extract)
{
	if (INFO_to_extract.size() == 0)
		return;

	printLOG("Outputting INFO for each site\n");
	unsigned int ui;
	string output = output_file_prefix + ".INFO";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open INFO Output File: " + output, 3);

	out << "CHROM\tPOS";
	for (ui=0; ui<INFO_to_extract.size(); ui++)
		out << "\t" << INFO_to_extract[ui];
	out << endl;

	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(false, false, true);

		out << e.get_CHROM() << "\t" << e.get_POS();

		for (ui=0; ui<INFO_to_extract.size(); ui++)
		{
			out << "\t" << e.get_INFO_value(INFO_to_extract[ui]);
		}
		out << endl;
	}

	out.close();
}

void vcf_file::output_genotype_depth(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output Genotype Depth Statistics.");

	printLOG("Outputting Depth for Each Genotype\n");
	string output = output_file_prefix + ".gdepth";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open Genotype Depth Output File: " + output, 7);

	out << "CHROM\tPOS";
	unsigned int ui;
	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		out << "\t" << indv[ui];
	}
	out << endl;

	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry();

		out << e.get_CHROM() << "\t" << e.get_POS();

		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.parse_genotype_entry(ui, false, false, true);
			out << "\t" << e.get_indv_DEPTH(ui);
		}
		out << endl;
	}
	out.close();
}

void vcf_file::output_TsTv(const string &output_file_prefix, int bin_size)
{
	printLOG("Outputting Ts/Tv Information\n");

	unsigned int N_alleles;

	map<string, unsigned int> model_to_idx;
	model_to_idx["AC"] = 0;
	model_to_idx["AG"] = 1;
	model_to_idx["AT"] = 2;
	model_to_idx["CG"] = 3;
	model_to_idx["CT"] = 4;
	model_to_idx["GT"] = 5;

	unsigned int s;
	map<string, int> max_pos;
	string vcf_line;
	string CHROM;
	vcf_entry e(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_basic_entry();

			CHROM = e.get_CHROM();

			if (max_pos.find(CHROM) != max_pos.end())
			{
				if (e.get_POS() > max_pos[CHROM])
					max_pos[CHROM] = e.get_POS();
			}
			else
				max_pos[CHROM] = e.get_POS();
		}
	}

	map<string, int>::iterator it;

	unsigned int N_bins;
	map<string, vector<int> > Ts_counts;
	map<string, vector<int> > Tv_counts;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		N_bins = (unsigned int)((max_pos[CHROM] + bin_size) / double(bin_size));
		Ts_counts[CHROM].resize(N_bins, 0);
		Tv_counts[CHROM].resize(N_bins, 0);
	}


	vector<unsigned int> model_counts(6,0);
	double C = 1.0 / double(bin_size);
	unsigned int idx;

	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		N_alleles = e.get_N_alleles();

		if (N_alleles != 2)
			continue;

		string model = e.get_REF() + e.get_ALT_allele(0);
		sort(model.begin(), model.end());

		CHROM = e.get_CHROM();
		idx = (unsigned int)(e.get_POS() * C);

		if (model_to_idx.find(model) != model_to_idx.end())
		{
			model_counts[model_to_idx[model]]++;
			switch (model_to_idx[model])
			{
			case 1:
			case 4:
				Ts_counts[CHROM][idx]++;
				break;
			case 0:
			case 2:
			case 3:
			case 5:
				Tv_counts[CHROM][idx]++;
				break;
			default:
				error("Unknown idx\n");
			}
		}
		else
		{
			error("Unknown model type\n");
		}
	}

	string output = output_file_prefix + ".TsTv";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open TsTv Output File: " + output, 7);

	out << "CHROM\tBinStart\tSNP_count\tTs/Tv" << endl;
	double ratio;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		for (s=0; s<Ts_counts[CHROM].size(); s++)
		{
			ratio = 0.0;
			if (Tv_counts[CHROM][s] != 0)
				ratio = double(Ts_counts[CHROM][s]) / Tv_counts[CHROM][s];
			out << CHROM << "\t" << s*bin_size << "\t" << Ts_counts[CHROM][s]+Tv_counts[CHROM][s] << "\t" << ratio << endl;
		}
	}


	out.close();

	output = output_file_prefix + ".TsTv.summary";
	out.open(output.c_str());
	if (!out.is_open())
		error("Could not open TsTv Summary Output File: " + output, 7);

	out << "MODEL\tCOUNT" << endl;
	out << "AC\t" << model_counts[0] << endl;
	out << "AG\t" << model_counts[1] << endl;
	out << "AT\t" << model_counts[2] << endl;
	out << "CG\t" << model_counts[3] << endl;
	out << "CT\t" << model_counts[4] << endl;
	out << "GT\t" << model_counts[5] << endl;
	unsigned int Ts = model_counts[1] + model_counts[4];
	unsigned int Tv = model_counts[0] + model_counts[2] + model_counts[3] + model_counts[5];
	out << "Ts\t" << Ts << endl;
	out << "Tv\t" << Tv << endl;

	printLOG("Ts/Tv ratio: " + dbl2str(double(Ts)/Tv, 3) + "\n");

	out.close();
}

void vcf_file::output_site_depth(const string &output_file_prefix, bool output_mean)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output Site Depth Statistics.");

	printLOG("Outputting Depth for Each Site\n");
	string output = output_file_prefix + ".ldepth";
	if (output_mean)
		output += ".mean";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open Site Depth Output File: " + output, 7);

	out << "CHROM\tPOS\t";
	if (output_mean)
		out << "MEAN_DEPTH" << endl;
	else
		out << "SUM_DEPTH" << endl;
	unsigned int ui;

	unsigned int N_include_indv = 0;
	if (output_mean)
		N_include_indv = N_kept_individuals();
	else
		N_include_indv = 1;

	int depth;
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry();

		out << e.get_CHROM() << "\t" << e.get_POS() << "\t";

		unsigned int sum=0;
		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.parse_genotype_entry(ui, false, false, true);
			depth = e.get_indv_DEPTH(ui);
			if (depth >= 0)
				sum += depth;
		}

		if (output_mean)
			out << double(sum) / N_include_indv << endl;
		else
			out << sum << endl;
	}
	out.close();
}


// Calculate, and output, Fst using the formula outlined in HapMap I
// Namely:
// Fst = 1 - (Pi_within / Pi_combined)
// where
// Pi_within = sum_j(nchoosek(n_j,2) * sum_i(2*n_ij * x_ij * (1-x_ij) / (n_ij -1))) / sum_j(nchoosek(n_j,2))
// and
// Pi_between = sum_i(2*n_i*x_i*(1-x_i) / (n_i - 1))
// where j is the population index, and i is the SNP index
void vcf_file::output_fst(const string &output_file_prefix, vcf_file &vcf_fst)
{
	printLOG("Outputting Fst estimates (for bi-allelic only)\n");

	string output = output_file_prefix + ".fst";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open Fst Output File: " + output, 7);

	out << "CHROM\tPOS\tFST" << endl;

	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;

	return_site_union(vcf_fst, CHROMPOS_to_filepos_pair);

	string vcf_line;

	int n_1, n_2, n_1_choose_2 = 0, n_2_choose_2=0;
	int last_n_1=-1, last_n_2=-1;

	unsigned int n_i1, n_i2, n_iT;
	int N_alleles1, N_alleles2;
	vector<int> allele_counts1, allele_counts2;
	double x_i1, x_i2, x_iT;
	int POS;
	int s1, s2;

	double tmp1, tmp2, tmpT;
	double sum1=0.0, sum2=0.0, sumT=0.0;
	double Fst;
	string CHROM;

	unsigned int N_intersecting_sites = 0;
	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		if ((s1 == -1) || (s2 == -1))
			continue;

		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		get_vcf_entry(s1, vcf_line);
		vcf_entry e1(N_indv, vcf_line);
		vcf_fst.get_vcf_entry(s2, vcf_line);
		vcf_entry e2(vcf_fst.N_indv, vcf_line);

		e1.parse_basic_entry(true);
		e2.parse_basic_entry(true);

		// Check sites have same alternative alleles
		N_alleles1 = e1.get_N_alleles();
		N_alleles2 = e2.get_N_alleles();

		if ((N_alleles1 > 2) || (N_alleles2 > 2))
			continue;

		if ((N_alleles1 == 1) && (N_alleles2 == 1))
			continue;

		if ((N_alleles1 == 2) && (N_alleles2 == 2))
			if (e1.get_ALT_allele(0) != e2.get_ALT_allele(0))
				continue;

		e1.parse_genotype_entries(true);
		e2.parse_genotype_entries(true);

		// Calculate allele frequencies
		e1.get_allele_counts(allele_counts1, n_i1, include_indv, include_genotype[s1]);
		e2.get_allele_counts(allele_counts2, n_i2, vcf_fst.include_indv, vcf_fst.include_genotype[s2]);

		if ((n_i1 == 0) || (n_i2 == 0))
			continue;

		n_1 = e1.get_N_chr(include_indv, include_genotype[s1]);
		n_2 = e2.get_N_chr(include_indv, include_genotype[s2]);

		if (last_n_1 != -1)
		{
			if ((n_1 != last_n_1) || (n_2 != last_n_2))
				error("Cannot mix sites with different ploidy. Are you including sex-chromosomes?\n");
		}
		else
		{
			last_n_1 = n_1;
			last_n_2 = n_2;
		}

		n_1_choose_2 = n_1 * (n_1 - 1) / 2;
		n_2_choose_2 = n_2 * (n_2 - 1) / 2;

		N_intersecting_sites++;

		x_i1 = allele_counts1[0] / double(n_i1);
		x_i2 = allele_counts2[0] / double(n_i2);
		n_iT = (n_i1 + n_i2);
		x_iT = (allele_counts1[0] + allele_counts2[0]) / double(n_iT);

		tmp1 = 2 * (n_i1 / (n_i1 - 1.0)) * x_i1 * (1-x_i1);
		tmp2 = 2 * (n_i2 / (n_i2 - 1.0)) * x_i2 * (1-x_i2);
		tmpT = 2 * (n_iT / (n_iT - 1.0)) * x_iT * (1-x_iT);

		Fst = 1.0 - (((n_1_choose_2 * tmp1) + (n_2_choose_2 * tmp2)) / (n_1_choose_2 + n_2_choose_2) / tmpT);

		out << CHROM << "\t" << POS << "\t" << Fst << endl;

		sum1 += tmp1;
		sum2 += tmp2;
		sumT += tmpT;

		last_n_1 = n_1; last_n_2 = n_2;
	}

	Fst = 1.0 - (((n_1_choose_2 * sum1) + (n_2_choose_2 * sum2)) / (n_1_choose_2 + n_2_choose_2) / sumT);

	printLOG("Found " + int2str(N_intersecting_sites) + " intersecting sites\n");
	printLOG("Fst = " + dbl2str(Fst, 6) + "\n");

	out.close();
}

void vcf_file::output_FORMAT_information(const string &output_file_prefix, const string &FORMAT_id)
{
	if (FORMAT_id == "")
		return;

	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output FORMAT information.");

	printLOG("Outputting FORMAT information for " + FORMAT_id + "\n");
	string output = output_file_prefix + "." + FORMAT_id + ".FORMAT";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open FORMAT Output File: " + output, 7);

	out << "CHROM\tPOS";
	unsigned int ui;
	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == true)
			out << "\t" << indv[ui];
	}
	out << endl;

	string vcf_line;
	vcf_entry e(N_indv);
	string FORMAT_out;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry();
		e.parse_full_entry(true);

		if (e.FORMAT_id_exists(FORMAT_id) == false)
			continue;

		out << e.get_CHROM() << "\t" << e.get_POS();

		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.read_indv_generic_entry(ui, FORMAT_id, FORMAT_out);
			out << "\t" << FORMAT_out;
		}
		out << endl;
	}


	out.close();
}

