__all__ = [ "Cmd" ]

import subprocess
import image
import tempfile

class Cmd(object):
    def __init__( self, program, *args, **kwargs ):
        self.program = program
        self.arguments = list(args) # list
        self.parameters = kwargs # dict
    def _ImageToFile( self, img ):
        (filehandle, tmp_file) = tempfile.mkstemp(suffix=".nii",dir=image.tmp_dir)
        image.imwrite( tmp_file, img )
        return tmp_file
    def _to_list(self):
        cmd_list = [ self.program ]
        if self.arguments is not None:
            for arg in self.arguments:
                if isinstance( arg, image.Image):
                    arg = self._ImageToFile(arg)
                cmd_list.append( arg )
        if self.parameters is not None:
            for name, value in self.parameters.items():
                cmd_list.append( '-' + name )
                if isinstance(value, str):
                    cmd_list.append(value)
                else:
                    for arg in value:
                        if isinstance( arg, image.Image):
                            arg = self._ImageToFile(arg)
                        cmd_list.append( arg )
        return cmd_list
    def __str__( self ):
        return ' '.join(self._to_list())
    def __repr__( self ):
        return str(self)     
    def __call__( self,  *args, **kwargs ):
        if args is not None:
            self.arguments.extend( args )
        if kwargs is not None:
            self.parameters.update( kwargs )
        proc = subprocess.Popen( self._to_list() )
        (out, err) = proc.communicate()
        if out is not None:
            print out
        if err is not None:
            print err
        return
        
