import importlib
import os
import tempfile
import traceback

import Pyro5.api

from df.data_transfer import DataFunctionRequest, DataFunction


@Pyro5.api.expose
class DataFunctions(object):

    def __init__(self):
        pass

    def data_function(self, request_json: str) -> str:
        try:
            request: DataFunctionRequest = DataFunctionRequest.parse_raw(request_json)
            id = request.id

            print(f'Processing request {id}')
            if request.script:
                exec(request.script, globals())
                response = execute(request)
            else:
                class_name = request.serviceName
                module = importlib.import_module(f'df.{class_name}')
                class_ = getattr(module, class_name)
                df: DataFunction = class_()
                response = df.execute(request)
            response_json = response.json()
            return response_json
        except Exception as ex:
            traceback.print_exc()
            print(f'Exception handling request {request.id}')
            raise ex


if __name__ == '__main__':
    if 'ENTREZ_EMAIL' not in os.environ:
        os.environ['ENTREZ_EMAIL'] = 'info@glysade.com'
    wd = os.path.join(tempfile.gettempdir(), 'pyro_df')
    if not os.path.exists(wd):
        os.mkdir(wd)
    os.chdir(wd)
    print(f'Sever working directory {wd}')
    # daemon = Pyro5.api.Daemon(host='0.0.0.0', port=9101)
    daemon = Pyro5.api.Daemon(host='127.0.0.1', port=9101)
    uri = daemon.register(DataFunctions, objectId='glysade.pyserver')
    print(f'Server ready URI {uri}')
    daemon.requestLoop()
