/*
Author: Brian C. Dandurand (2017-2018)
Based on and modified from code developed by Jeffrey Christiansen, Brian Dandurand, and Fabricio Oliveira
at RMIT in Melbourne Australia with funding under project ARC DP 140100985 during 2015-2017.
CIs of that projects were Prof. Andrew Eberhard, Prof. Natashia Boland, and PI Prof. Jeffrey Linderoth.
*/

#ifndef PSCGSCEN_H
#define PSCGSCEN_H

#include "OsiSolverInterface.hpp"
#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"


class PSCGScen{
public:
    PSCGScen(SmiScnModel &smi, OsiSolverInterface &smi_osi, OsiSolverInterface &sp_osi, int scn):smi(smi),smi_osi(smi_osi),sp_osi(sp_osi),scn(scn){
        SmiScnNode *node = smi.getLeafNode(scn);
	n_stg = 0;
        while (node != NULL){
            cols_by_stg.push_back(vector<int>());
            orig_col_ubs_by_stg.push_back(vector<double>());
            orig_col_lbs_by_stg.push_back(vector<double>());
            orig_cfts_by_stg.push_back(vector<double>());

            rows_by_stg.push_back(vector<int>());
            orig_row_ubs_by_stg.push_back(vector<double>());
            orig_row_lbs_by_stg.push_back(vector<double>());

            stg_probs.push_back(node->getProb());

            node = node->getParent();
            n_stg++;
        }
        int stg = n_stg;
        node = smi.getLeafNode(scn);
        while (node != NULL){
            stg--;
            for(int jj=node->getColStart(); jj<node->getColStart()+node->getNumCols(); ++jj){
                cols_by_stg[stg].push_back(jj); 
            	orig_cfts_by_stg[stg].push_back( (smi_osi.getObjCoefficients())[jj] );
                orig_col_lbs_by_stg[stg].push_back( (smi_osi.getColLower())[jj] );
                orig_col_ubs_by_stg[stg].push_back( (smi_osi.getColUpper())[jj] );
            }
            for(int jj=node->getRowStart(); jj<node->getRowStart()+node->getNumRows(); ++jj){
                rows_by_stg[stg].push_back(jj); 
                orig_row_lbs_by_stg[stg].push_back( (smi_osi.getRowLower())[jj] );
                orig_row_ubs_by_stg[stg].push_back( (smi_osi.getRowUpper())[jj] );
            }
            node = node->getParent();
        }
        for( int stg=0; stg<n_stg; stg++){
            for( int row=0; row<rows_by_stg[stg].size(); row++){
                all_scn_rows.push_back(rows_by_stg[stg][row]);
                all_orig_row_ubs.push_back(orig_row_ubs_by_stg[stg][row]);
                all_orig_row_lbs.push_back(orig_row_lbs_by_stg[stg][row]);
            }
            for( int col=0; col<cols_by_stg[stg].size(); col++){
                all_scn_cols.push_back(cols_by_stg[stg][col]);

                all_orig_cfts.push_back(orig_cfts_by_stg[stg][col]);
                all_orig_col_ubs.push_back(orig_col_ubs_by_stg[stg][col]);
                all_orig_col_lbs.push_back(orig_col_lbs_by_stg[stg][col]);

                all_scn_col_probs.push_back(stg_probs[stg]);
            }
        }
        nrows = all_scn_rows.size();
        ncols = all_scn_cols.size();

	for(int col=0; col<ncols; col++){
            if(smi_osi.isInteger(  all_scn_cols[col]  )){
                int_indx.push_back(col);
            }
        }

        constructScnSubmatrix( ); //extract the scenario-specific submatrix of smi_osi into sp_osi
        setupOsi();
    }
    ~PSCGScen(){
    }
    int getNumStages(){return cols_by_stg.size();}
    int getNumColsByStage(int stg){return cols_by_stg[stg].size();}
    int getNumCols(){return ncols;}
    int getNumRowsByStage(int stg){return rows_by_stg[stg].size();}
    int getNumRows(){return nrows;}
    vector<int> getIntInds(){return int_indx;}
    void enforceIntegrality(){
        for (unsigned int nn = 0; nn < int_indx.size(); nn++) {
            sp_osi.setInteger(int_indx[nn]);
        }
    }
    void modifyObjCft(int col, double val){
            sp_osi.setObjCoeff(col, val );
    }
    void resetObjCfts(){
        for( int col=0; col < ncols; col++){
            sp_osi.setObjCoeff(col,  all_scn_col_probs[col]*all_orig_cfts[col]  );
        }
    }
    void resetColBds(){
        for( int col=0; col < ncols; col++){
            sp_osi.setColLower(col, all_orig_col_lbs[col] );
            sp_osi.setColUpper(col, all_orig_col_ubs[col] ); 
        }
    }
    void resetRowBds(){
        for( int row=0; row < nrows; row++){
            sp_osi.setRowLower(row, all_orig_row_lbs[row] );
            sp_osi.setRowUpper(row, all_orig_row_ubs[row] ); 
        }
    }

protected:

    int scn;
    int n_stg;
    int ncols;
    int nrows;
    SmiScnModel &smi;
    OsiSolverInterface &smi_osi;
    OsiSolverInterface &sp_osi;
    vector<vector<int> > cols_by_stg;
    vector<vector<int> > rows_by_stg;
    vector<int> all_scn_rows;
    vector<int> all_scn_cols;
    vector<double> stg_probs;
    vector<double> all_scn_col_probs;
    vector<int> int_indx;
    vector< vector<int> > int_indx_by_stg;

    vector<double> all_orig_cfts;
    vector<vector<double> > orig_cfts_by_stg;

    vector<double> all_orig_col_ubs;
    vector<vector<double> > orig_col_ubs_by_stg;
    vector<double> all_orig_col_lbs;
    vector<vector<double> > orig_col_lbs_by_stg;
    vector<double> all_orig_row_ubs;
    vector<vector<double> > orig_row_ubs_by_stg;
    vector<double> all_orig_row_lbs;
    vector<vector<double> > orig_row_lbs_by_stg;

    CoinPackedMatrix sp_mat;

private:
    //The following methods are helper methods for the constructor and are not meant to be part of the public interface.
    void constructScnSubmatrix( ){
        const CoinPackedMatrix *smi_matptr = smi_osi.getMatrixByCol();
        sp_mat.clear();
        sp_mat.setDimensions(nrows,ncols);

    	int from_row,from_col;
        for( int row=0; row < nrows; row++ ){
            for( int col=0; col < ncols; col++){
	        from_col = all_scn_cols[col];
	        from_row = all_scn_rows[row];
	        double elmt = smi_matptr->getCoefficient(from_row,from_col);
	        sp_mat.modifyCoefficient(row,col,elmt);
            }
        }
        return;
    }
    void setupOsi(){
        sp_osi.loadProblem(sp_mat,NULL,NULL,NULL,NULL,NULL);
	resetColBds();
	resetRowBds();
	resetObjCfts();
        sp_osi.setObjSense(smi_osi.getObjSense());  // Set objective sense, MIN=1.0
    }

};

//TODO: Extend the above class as an abstract class with certain inherited interfaces corresponding to the specific Osi 

#endif
