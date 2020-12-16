// Hacky code to access protected methods
// This converts protected methods/vars into public ones, so that they can be
// used in PyROOT (which can't see protected things)
// 
// IMPORTANT: don't start with a block comment...wont work

#ifndef My_TUnfoldDensity
#define My_TUnfoldDensity

// #include "TUnfoldDensity.h"

class MyTUnfoldDensity : public TUnfoldDensity {
public:
    MyTUnfoldDensity(void);   // constructor for derived classes, do nothing

    MyTUnfoldDensity(const TH2 *hist_A,
                     EHistMap histmap,
                     ERegMode regmode = kRegModeCurvature,
                     EConstraint constraint=kEConstraintArea,
                     EDensityMode densityMode=kDensityModeBinWidthAndUser,
                     const TUnfoldBinning *outputBins=0,
                     const TUnfoldBinning *inputBins=0,
                     const char *regularisationDistribution=0,
                     const char *regularisationAxisSteering="*[UOB]"
    ) {
        
        TUnfoldDensity(hist_A,
                       histmap,
                       regmode,
                       constraint,
                       densityMode,
                       outputBins,
                       inputBins,
                       regularisationDistribution,
                       regularisationAxisSteering);
       std::cout << "Setup (My)TUnfoldDensity" << std::endl;
    }

    // Here are methods that were protected in TUnfoldDensity, that are now public for us to use
 
    virtual TString GetOutputBinName(Int_t iBinX) const { return TUnfold::GetOutputBinName(iBinX); } // name a bin
    // multiply sparse and non-sparse matrix
    TMatrixDSparse *MultiplyMSparseM(const TMatrixDSparse *a,const TMatrixD *b) const { return TUnfold::MultiplyMSparseM(a, b); } 
    // multiply sparse and sparse matrix
    TMatrixDSparse *MultiplyMSparseMSparse(const TMatrixDSparse *a,const TMatrixDSparse *b) const { return TUnfold::MultiplyMSparseMSparse(a, b); } 
    // multiply transposed sparse and sparse matrix
    TMatrixDSparse *MultiplyMSparseTranspMSparse(const TMatrixDSparse *a,const TMatrixDSparse *b) const { return TUnfold::MultiplyMSparseTranspMSparse(a, b); } 
    // calculate M_ij = sum_k [m1_ik*m2_jk*v[k] ]. the pointer v may be zero (means no scaling).
    TMatrixDSparse *MultiplyMSparseMSparseTranspVector (const TMatrixDSparse *m1,const TMatrixDSparse *m2, const TMatrixTBase<Double_t> *v) const { return TUnfold::MultiplyMSparseMSparseTranspVector(m1, m2, v); } 
    // invert symmetric (semi-)positive sparse matrix
    TMatrixDSparse *InvertMSparseSymmPos(const TMatrixDSparse *A, bool doPseudoInverse) const { 
        if (doPseudoInverse) { 
            Int_t rank = 1; 
            return TUnfold::InvertMSparseSymmPos(A, &rank); 
        } else { 
            return TUnfold::InvertMSparseSymmPos(A, nullptr); 
        } 
    }
    // replacement for dest += f*src
    void AddMSparse(TMatrixDSparse *dest,Double_t f,const TMatrixDSparse *src) const { return TUnfold::AddMSparse(dest, f, src); } 
    // create a TMatrixDSparse from an array
    TMatrixDSparse *CreateSparseMatrix(Int_t nrow,Int_t ncol,Int_t nele,Int_t *row,Int_t *col,Double_t *data) const { return TUnfold::CreateSparseMatrix(nrow, ncol, nele, row, col, data); } 
    TMatrixDSparse *CreateSparseMatrix(Int_t nrow,Int_t ncol,Int_t nele,std::vector<Int_t> row,std::vector<Int_t> col,std::vector<Double_t> data) const { return TUnfold::CreateSparseMatrix(nrow, ncol, nele, &row[0], &col[0], &data[0]); } 

    /// returns internal number of output (truth) matrix rows
    inline Int_t GetNx(void) const { return TUnfold::GetNx(); }
    /// converts truth histogram bin number to matrix row
    inline Int_t GetRowFromBin(int ix) const { return TUnfold::GetRowFromBin(ix); }
    /// converts matrix row to truth histogram bin number
    inline Int_t GetBinFromRow(int ix) const { return TUnfold::GetBinFromRow(ix); }
    /// returns the number of measurement bins
    inline Int_t GetNy(void) const { return TUnfold::GetNy(); }
    /// vector of the unfolding result
    inline const TMatrixD *GetX(void) const { return TUnfoldDensity::GetX(); }
    /// covariance matrix of the result
    inline const TMatrixDSparse *GetVxx(void) const { return TUnfoldDensity::GetVxx(); }
    /// inverse of covariance matrix of the result
    inline const TMatrixDSparse *GetVxxInv(void) const { return TUnfoldDensity::GetVxxInv(); }
    /// vector of folded-back result
    inline const TMatrixDSparse *GetAx(void) const { return TUnfoldDensity::GetAx(); }
    /// matrix of derivatives dx/dy
    inline const TMatrixDSparse *GetDXDY(void) const { return TUnfoldDensity::GetDXDY(); }
    /// matrix contributions of the derivative dx/dA
    inline const TMatrixDSparse *GetDXDAM(int i) const { return TUnfoldDensity::GetDXDAM(i); }
    /// vector contributions of the derivative dx/dA
    inline const TMatrixDSparse *GetDXDAZ(int i) const { return TUnfoldDensity::GetDXDAZ(i); }
    /// matrix E<sup>-1</sup>, using internal bin counting
    inline const TMatrixDSparse *GetEinv(void) const { return TUnfoldDensity::GetEinv(); }
    /// matrix E, using internal bin counting
    inline const TMatrixDSparse *GetE(void) const { return TUnfoldDensity::GetE(); }
    /// vector of the input data y (background-subtracted)
    inline const TMatrixD *GetY(void) const { return fY; }
    /// covariance matrix of the data y
    inline const TMatrixDSparse *GetVyy(void) const { return fVyy; }
    /// inverse of covariance matrix of the data y
    inline const TMatrixDSparse *GetVyyInv(void) const { return TUnfoldDensity::GetVyyInv(); }
    /// response matrix A
    inline const TMatrixDSparse *GetA(void) const { return fA; }
    /// mapping array to hist
    inline const TArrayI & GetXToHist(void) const { return fXToHist; }
    /// mapping hist to array
    inline const TArrayI & GetHistToX(void) const { return fHistToX; }
    
   // return an error matrix as histogram
    void ErrorMatrixToHist(TH2 *ematrix, const TMatrixDSparse *emat) const { TUnfold::ErrorMatrixToHist(ematrix, emat, nullptr, true); }

    /// delete matrix and invalidate pointer
    static void DeleteMatrix(TMatrixD **m) { TUnfold::DeleteMatrix(m); }
    /// delete sparse matrix and invalidate pointer
    static void DeleteMatrix(TMatrixDSparse **m) { TUnfold::DeleteMatrix(m); }

    /// determine total error matrix on the vector Ax
    TMatrixDSparse *GetSummedErrorMatrixYY(void) { return TUnfoldSys::GetSummedErrorMatrixYY(); }
    
    /// determine total error matrix on the vector x
    TMatrixDSparse *GetSummedErrorMatrixXX(void) { return TUnfoldSys::GetSummedErrorMatrixXX(); }

    Bool_t AddRegularisationCondition(Int_t i0, Double_t f0, Int_t i1=-1, Double_t f1=0., Int_t i2=-1, Double_t f2=0.) { return TUnfold::AddRegularisationCondition(i0,f0,i1,f1,i2,f2); } // add regularisation condition for a triplet of output bins
    Bool_t AddRegularisationCondition(Int_t nEle, const Int_t *indices, const Double_t *rowData) { return TUnfold::AddRegularisationCondition(nEle, indices, rowData); } // add a regularisation condition

    // FIXME: this shouldn't need an arg?
    vector<Int_t> CheckInput(const TH1 *hist_y=nullptr)
    {
        vector<Int_t> ignoredInputBins = {};

        Double_t oneOverZeroError = 0.0;

        Int_t *rowVyy1=new Int_t[GetNy()];
        Int_t *colVyy1=new Int_t[GetNy()];
        Double_t *dataVyy1=new Double_t[GetNy()];
        Double_t *dataVyyDiag=new Double_t[GetNy()];

        Int_t nVarianceZero=0;
        Int_t nVyy1=0;
        if (hist_y == nullptr && fVyy == nullptr) {
            throw std::runtime_error("hist_y and fVyy are empty, either specify hist_y or call SetInput() first");
        }
        for (Int_t iy = 0; iy < GetNy(); iy++) {
            // diagonals
            Double_t dy;
            if (hist_y != nullptr)
                dy = hist_y->GetBinError(iy + 1);
            else
                dy = (*fVyy)(iy, iy);
            Double_t dy2=dy*dy;
            if (dy2 <= 0.0) {
                nVarianceZero++;
            }
            rowVyy1[nVyy1] = iy;
            colVyy1[nVyy1] = 0;
            dataVyyDiag[iy] = dy2;
            if(dy2>0.0) {
                dataVyy1[nVyy1++] = dy2;
            }
        }

        TMatrixDSparse *vecV=CreateSparseMatrix(GetNy(),1,nVyy1,rowVyy1,colVyy1, dataVyy1);

        delete[] rowVyy1;
        delete[] colVyy1;
        delete[] dataVyy1;

        // simple check whether unfolding is possible, given the matrices fA and fV
        TMatrixDSparse *mAtV=MultiplyMSparseTranspMSparse(fA,vecV);
        DeleteMatrix(&vecV);
        Int_t nError2=0;
        //mAtV is a nx1 matrix
        for (Int_t i = 0; i<mAtV->GetNrows();i++) {
            if(mAtV->GetRowIndexArray()[i]==mAtV->GetRowIndexArray()[i+1]) {
                nError2 ++;
            }
        }
        if(nVarianceZero) {
            if(nVarianceZero>1) {
                Warning("CheckInput","%d/%d input bins have zero error,"
                    " and are ignored.",nVarianceZero,GetNy());
            } else {
                Warning("CheckInput","One input bin has zero error,"
                    " and is ignored.");
            }
        }
        // printout mAtV vector
        for (Int_t irow = 0; irow < mAtV->GetNrows(); irow++) {
          const Int_t sIndex = mAtV->GetRowIndexArray()[irow];
          const Int_t eIndex = mAtV->GetRowIndexArray()[irow+1];
          if (sIndex == eIndex) {
            cout << "sIndex==eIndex for sIndex= " << sIndex << endl;
          }
          for (Int_t index = sIndex; index < eIndex; index++) {
            const Int_t icol = mAtV->GetColIndexArray()[index];
            const Double_t data = mAtV->GetMatrixArray()[index];
            // printf("mAtV(%d,%d) = %.4e\n",irow,icol,data);
          }
        }

        if(nError2>0) {
            // check whether data points with zero error are responsible
            if(oneOverZeroError<=0.0) {
                // a_rows[i] has the start index (in a_cols, a_data) 
                // for row i, a_rows[i+1] is the end index + 1
                // a_cols[a_rows[i]]...a_cols[a_rows[i+1]-1] 
                // is the set of column indices that have !=0 data
                // a_data[a_rows[i]]...a_data[a_rows[i+1]-1] 
                // is the set of corresponding data
                cout << "fA NRows: " << fA->GetNrows() << endl;
                TMatrixDSparse *aT = new TMatrixDSparse(fA->GetNcols(), fA->GetNrows());
                aT->Transpose(*fA);
                cout << "aT NRows: " << aT->GetNrows() << endl;

                const Int_t *a_rows = aT->GetRowIndexArray(); 
                const Int_t *a_cols = aT->GetColIndexArray();
                const Double_t *a_data = aT->GetMatrixArray();

                // for each row in the resultant mAtV, look for 0 entries.
                // for each of those, look at all the A^T and dataVyyDiag pairs
                // (i.e. iterate over columns) and see which are 0
                for (Int_t row = 0; row <mAtV->GetNrows();row++) {
                    // printout every row's sum entries
                    // printf("++++ ROW %d\n", row);
                    // const Int_t sIndex = a_rows[row];
                    // const Int_t eIndex = a_rows[row+1];
                    // for (Int_t i=sIndex; i<eIndex; i++) {
                    //     Int_t col = a_cols[i];
                    //     cout << "dataVyyDiag[" << col << "]: " << dataVyyDiag[col] << endl;
                    //     printf("a^T(%d,%d) = %.4e\n", row, col, a_data[i]);
                    // }

                    if (mAtV->GetRowIndexArray()[row] == mAtV->GetRowIndexArray()[row+1]) {
                        TString binlist("no data to constrain output bin ");
                        binlist += GetOutputBinName(fXToHist[row]);
                        binlist +=" depends on ignored input bins ";

                        const Int_t sIndex = a_rows[row];
                        const Int_t eIndex = a_rows[row+1];
                        for (Int_t i=sIndex; i<eIndex; i++) {
                            Int_t col = a_cols[i];
                            if (dataVyyDiag[col] > 0.) continue;
                            // cout << "dataVyyDiag[" << col << "]: " << dataVyyDiag[col] << endl;
                            // printf("a^T(%d,%d) = %.4e\n", row, col, a_data[i]);
                            binlist +=" ";
                            binlist += col+1; // +1 for ROOT binning
                            ignoredInputBins.push_back(col+1);
                        }
                        // inhereted from TObject
                        Warning("CheckInput",binlist);
                    }
                }
            }
            
            // TUnfold original version of the same thing that I can't understand :(
            // just to check I did it right
            if(oneOverZeroError<=0.0) {
                const Int_t *a_rows=fA->GetRowIndexArray();
                const Int_t *a_cols=fA->GetColIndexArray();
                for (Int_t col = 0; col <mAtV->GetNrows();col++) {
                    if(mAtV->GetRowIndexArray()[col]==mAtV->GetRowIndexArray()[col+1]) {
                        TString binlist("no data to constrain output bin ");
                        binlist += GetOutputBinName(fXToHist[col]);
                        binlist +=" depends on ignored input bins ";
                        for(Int_t row=0;row<fA->GetNrows();row++) {
                            if(dataVyyDiag[row]>0.0) continue;
                            for(Int_t i=a_rows[row];i<a_rows[row+1];i++) {
                                if(a_cols[i]!=col) continue;
                                binlist +=" ";
                                binlist +=row+1;
                            }
                        }  
                        Warning("SetInput", binlist);
                    }
                }
            }

            if(nError2>1) {
                Error("CheckInput","%d/%d output bins are not constrained by any data.",
                    nError2,mAtV->GetNrows());
            } else {
                Error("CheckInput","One output bins is not constrained by any data.");
            }
        }
        DeleteMatrix(&mAtV);

        delete[] dataVyyDiag;

        return ignoredInputBins;
    }

    ClassDef(MyTUnfoldDensity, 1) //My Unfolding with density regularisation

};

ClassImp(MyTUnfoldDensity);

#endif