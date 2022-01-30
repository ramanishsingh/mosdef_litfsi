"""Configuration of the project enviroment. 
 
The environments defined in this module can be auto-detected. 
This helps to define environment specific behaviour in heterogenous 
environments. 
""" 
import flow 
 
__all__ = ['CoriEnvironment']
 
# class CoriEnvironment(flow.environment.SlurmEnvironment): 
#     hostname_pattern = 'master.cl.vanderbilt.edu' 
#     template = 'rahman.sh'
#     cores_per_node = 16
class RahmanEnvironment(flow.environment.DefaultPBSEnvironment):
    hostname_pattern = 'master.cl.vanderbilt.edu'
    template = 'rahman.vanderbilt.sh'
    cores_per_node = 16
    
    @classmethod
    def add_args(cls, parser):
        super(flow.environment.DefaultPBSEnvironment, cls).add_args(parser)
        parser.add_argument(
               '-w', "--walltime", type=float, help="Walltime"
            )
