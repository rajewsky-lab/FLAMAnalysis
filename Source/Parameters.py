import yaml
import sys
from Source import Misc

"""
class Parameters()

- Load FLAMSeq Experiment and Software Parameters
- Define createdDirs
- Provide regDirs methods -> self.regDirs contains dict with keys and paths for dirs that have been created by
    pipeline session

"""

class Parameters():

    def __init__(self, yamlParameterPath):
        self.parameterPath = yamlParameterPath

        # Check if YAML file exists
        Misc.Misc.checkFile(self.parameterPath, 'Parameters YAML File')

        self.parameters = self.parseYaml(self.parameterPath)

        # regDirs contains all dirs that have been created / checked during pipeline session.
        self.regDirs = dict()

    def parseYaml(self, yamlParameterPath):
        with open(yamlParameterPath, 'r') as yamlStream:
            try:
                parameters = yaml.load(yamlStream)
            except yaml.YAMLError:
                sys.stderr.write("Error Reading Parameter YAML at: %s\n".format(yamlParameterPath))

        return parameters

    def getParameters(self):
        return self.parameters

    def getParametersByKey(self, key):
        try:
            return self.parameters[key]
        except KeyError:
            raise KeyError('No Data found for Key for Parameter')

    def registerDir(self, key, path):
    # Register dir path for key

        if key in self.regDirs:
            raise KeyError('Key for regDict already registered.')
        else:
            self.regDirs[key] = path

    def getRegDir(self, key):

        if not key in self.regDirs:
            return KeyError('Key not in regDict')
        else:
            return self.regDirs[key]

    def checkKey(self, key):
        if key in self.regDirs:
            return True
        else:
            return False

    def getExpName(self):
        try:
            return self.parameters['experimentName']
        except KeyError:
            raise KeyError("No Experiment Defined In parameters.yaml")