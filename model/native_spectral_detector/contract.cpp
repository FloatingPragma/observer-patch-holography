// Full-row contraction of an already compiled native local instrument.
// This performs no spectral, particle-number, amplitude or rank truncation.
#include <cstdint>
#include <cstring>

template<int D>
static void apply(double* state, const std::int64_t* bases, std::int64_t count,
                  const std::int64_t* offsets, int columns,
                  const double* matrix, double* density) {
    double old[D * 32 * 2];
    double rho[D * D * 2] = {};
    for (std::int64_t b = 0; b < count; ++b) {
        for (int r = 0; r < D; ++r)
            std::memcpy(old + 2*r*columns,
                        state + 2*(bases[b] | offsets[r])*columns,
                        2*columns*sizeof(double));
        for (int r = 0; r < D; ++r) {
            for (int s = r; s < D; ++s) {
                double real = 0., imag = 0.;
                for (int j = 0; j < columns; ++j) {
                    const double ar=old[2*(r*columns+j)], ai=old[2*(r*columns+j)+1];
                    const double br=old[2*(s*columns+j)], bi=old[2*(s*columns+j)+1];
                    real += ar*br+ai*bi;
                    imag += ai*br-ar*bi;
                }
                rho[2*(r*D+s)] += real;
                rho[2*(r*D+s)+1] += imag;
            }
            for (int j = 0; j < columns; ++j) {
                double real = 0., imag = 0.;
                for (int s = 0; s < D; ++s) {
                    const double ar=matrix[2*(r*D+s)], ai=matrix[2*(r*D+s)+1];
                    const double br=old[2*(s*columns+j)], bi=old[2*(s*columns+j)+1];
                    real += ar*br-ai*bi;
                    imag += ar*bi+ai*br;
                }
                const std::int64_t index=2*((bases[b] | offsets[r])*columns+j);
                state[index]=real;
                state[index+1]=imag;
            }
        }
    }
    for (int r=0;r<D;++r) for (int s=r;s<D;++s) {
        density[2*(r*D+s)]=rho[2*(r*D+s)];
        density[2*(r*D+s)+1]=rho[2*(r*D+s)+1];
        density[2*(s*D+r)]=rho[2*(r*D+s)];
        density[2*(s*D+r)+1]=-rho[2*(r*D+s)+1];
    }
}

extern "C" int oph_apply(double* state, const std::int64_t* bases, std::int64_t count,
                         const std::int64_t* offsets, int columns, int dimension,
                         const double* matrix, double* density) {
    if (columns<1 || columns>32 || count<1) return 1;
    switch (dimension) {
        case 2:apply<2>(state,bases,count,offsets,columns,matrix,density);break;
        case 4:apply<4>(state,bases,count,offsets,columns,matrix,density);break;
        case 8:apply<8>(state,bases,count,offsets,columns,matrix,density);break;
        default:return 2;
    }
    return 0;
}
