import importlib
import sys

from df.data_transfer import DataFunctionRequest

if __name__ == '__main__':
    in_file = sys.argv[1]
    out_file = sys.argv[2]

    with open(in_file, 'r') as fh:
        request_json = fh.read()

    request = DataFunctionRequest.parse_raw(request_json)  # type: DataFunctionRequest
    class_name = request.serviceName
    module = importlib.import_module(f'df.{class_name}')
    class_ = getattr(module, class_name)
    df = class_()  # type: DataFunction
    response = df.execute(request)
    response_json = response.json()

    with open(out_file, 'w') as fh:
        fh.write(response_json)
