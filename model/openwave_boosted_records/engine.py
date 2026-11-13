"""Standalone Taichi complete-sector propagation; no RER imports."""
import argparse
import json
from pathlib import Path
import numpy as np
import taichi as ti


@ti.data_oriented
class SectorKernel:
    def __init__(self,n,frame):
        self.n=n;self.frame=ti.Vector.field(2,ti.f64,shape=(96,48))
        self.frame.from_numpy(np.stack((frame.real,frame.imag),axis=-1))
        self.state=ti.Vector.field(2,ti.f64,shape=(n,48))
        self.energy=ti.field(ti.f64,shape=(n,48))
        self.native=ti.Vector.field(2,ti.f64,shape=(n,96))
        self.record=ti.Vector.field(2,ti.f64,shape=(n,96))

    @ti.kernel
    def evolve(self,time:ti.f64):
        for k,j in self.state:
            phase=time*self.energy[k,j];c=ti.cos(phase);s=ti.sin(phase)
            v=self.state[k,j];self.state[k,j]=ti.Vector([c*v[0]+s*v[1],c*v[1]-s*v[0]])

    @ti.kernel
    def probe(self):
        for k,a in self.native:
            out=ti.Vector([0.,0.])
            for j in range(48):
                u=self.frame[a,j];v=self.state[k,j]
                out+=ti.Vector([u[0]*v[0]-u[1]*v[1],u[0]*v[1]+u[1]*v[0]])
            self.native[k,a]=out
        for x,a in self.record:
            out=ti.Vector([0.,0.])
            for k in range(self.n):
                phase=2*np.pi*k*x/self.n;c=ti.cos(phase);s=ti.sin(phase);v=self.native[k,a]
                out+=ti.Vector([c*v[0]-s*v[1],c*v[1]+s*v[0]])/ti.sqrt(ti.cast(self.n,ti.f64))
            self.record[x,a]=out


def execute(packet):
    ti.init(arch=ti.cpu,default_fp=ti.f64,random_seed=0,cpu_max_num_threads=1)
    T=np.asarray(packet['native_positive_frame_real'])+1j*np.asarray(packet['native_positive_frame_imag'])
    masses=np.asarray(packet['masses_electron_units']);rows=[]
    for n in packet['sides']:
        length=8*np.pi/3;a=length/n;p=2*np.pi*np.fft.fftfreq(n,d=a)
        energy=np.sqrt(masses[None,:]**2+(2*np.sin(a*p/2)/a)[:,None]**2)
        cont=np.sqrt(masses[None,:]**2+p[:,None]**2)
        kernel=SectorKernel(n,T);kernel.energy.from_numpy(energy)
        for case in packet['cases']:
            initial=np.zeros((n,48),complex);j=case['mode'];weights=np.asarray(case['amplitudes_real'],float)
            weights/=np.linalg.norm(weights)
            for k,v in zip(case['momenta_in_grid_units'],weights):initial[k%n,j]=v
            kernel.state.from_numpy(np.stack((initial.real,initial.imag),axis=-1));kernel.evolve(case['time']);kernel.probe()
            raw=kernel.record.to_numpy();wave=raw[...,0]+1j*raw[...,1]
            independent=np.fft.ifft((initial*np.exp(-1j*case['time']*energy))@T.T,axis=0,norm='ortho')
            continuum=np.fft.ifft((initial*np.exp(-1j*case['time']*cont))@T.T,axis=0,norm='ortho')
            probabilities=.5*abs(wave)**2;reference=.5*abs(continuum)**2
            tv=float(np.sum(abs(probabilities-reference))/2)
            P=max(abs(p[k%n]) for k in case['momenta_in_grid_units']);bound=min(1,.5*case['time']*a*a*P**4/(24*masses[j]))
            floating=float(np.linalg.norm(wave-independent));assert floating<3e-12
            assert tv<=bound+3e-12 and abs(float(np.sum(probabilities))-.5)<3e-12
            reverse=np.fft.ifft((initial*np.exp(+1j*case['time']*energy))@T.T,axis=0,norm='ortho')
            mutation=float(np.sum(abs(probabilities-.5*abs(reverse)**2))/2)
            rows.append({'side':n,'case':case['name'],'complete_position_native_probe_probabilities':probabilities.tolist(),
                'preparation_failure_probability':.5,'Gauss_constraint_failure_probability':0.,
                'all_other_Fock_sector_probability':0.,'complete_probability_sum_diagnostic':float(np.sum(probabilities)+.5),
                'native_coordinate_count':96,'positive_source_mode_count':48,
                'independent_full_sector_state_error':floating,
                'same_input_continuum_complete_record_TV_diagnostic':tv,
                'continuum_record_TV_upper_diagnostic':bound,
                'reversed_time_control_complete_record_TV':mutation,
                'energy_mean_electron_units_diagnostic':float(np.sum(abs(initial)**2*energy)),
                'success_preparation_work_electron_units_diagnostic':.5*float(np.sum(abs(initial)**2*energy)),
                'clock_parameter_electron_inverse_mass':case['time']})
    assert min(r['reversed_time_control_complete_record_TV'] for r in rows if r['case']=='electron_packet')>.05
    import openwave
    return {'solver':'Taichi CPU float64 exact spectral phases and full native probe transform',
        'taichi_version':list(ti.__version__),'OpenWave_version':openwave.__version__,
        'rows':rows,'source_masses_all48':masses.tolist(),
        'all_native_probe_rows_and_zero_probability_complements_retained':True}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=execute(json.loads(args.source.read_text()));args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')


if __name__=='__main__':main()
