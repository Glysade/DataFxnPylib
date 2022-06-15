"""
================
ruse_resource.py
================

Copyright (C) 2017-2022 Glysade, LLC

Contains a class for manipulating ruse resources

"""


import os
from typing import Iterable, Dict, Optional

import requests
from ruse.util.frozen import Frozen
from ruse.util.util import ruse_server, RuseServerException, ruse_work_dir


class RuseResource(Frozen):
    """
    A class for managing a ruse resource.  Most of the attributes map to those stored in the ruse database

    Attributes:

        - owner: the name of the resource owner

        - resource_type: the ruse resource type

        - orig_filename: the original filename

        - name: resource name

        - meta: resource meta data

        - resource_id: resource identifier

        - binary_file: if True, resource is binary

        - owner_id: the id of the owner

        - resource_type_id: the id of the resource type

        - uuid: UUID of the resource

        - ext: extension of the resource

        - downloaded_file: path to local copy of downloaded resource

    """
    def __init__(self):
        self.owner = None  # type: str
        self.resource_type = None  # type: str
        self.orig_filename = None  # type: str
        self.name = None  # type: str
        self.meta = None  # type: str
        self.resource_id = None  # type: int
        self.binary_file = True

        self.owner_id = None  # type: int
        self.resource_type_id = None  # type: int
        self.uuid = None  # type: str
        self.ext = None  # type: str
        self.downloaded_file = None  # type: str

    @classmethod
    def _get_error(cls, response) -> str:
        """
        A class method that parses http response to get any json error

        :param response: http response
        :return: json error
        """

        error = "Unknown error"
        try:
            error = response.json()['error']
        finally:
            return error

    def upload_file(self, resource_type: str, orig_filename: str, owner: str = 'default', name: str = None,
                    meta: str = None, binary_file: bool = True):
        """
        Uploads a resource file to the Ruse server

        :param resource_type: the type of the resource
        :param orig_filename: a local file name for the resource
        :param owner: the owner of the resource (default "default")
        :param name: the name of the resource (default None)
        :param meta: meta data (default None)
        :param binary_file: True (default) if the file is binary
        :return: uploaded resource id
        """

        assert resource_type is not None
        assert orig_filename is not None

        self.owner = owner
        self.resource_type = resource_type
        self.orig_filename = orig_filename
        self.name = name
        self.meta = meta
        self.binary_file = binary_file

        if not self.name:
            self.name = resource_type

        # request a resource id
        fields = ['owner', 'resource_type', 'orig_filename', 'name', 'meta']
        resource_json = {attr: getattr(self, attr) for attr in fields if getattr(self, attr)}
        uri = "{}/resource".format(ruse_server())
        headers = {'Content-Type': 'application/json'}
        response = requests.post(uri, json=resource_json, headers=headers)
        if response.status_code != requests.codes.ok:
            raise RuseServerException(
                "upload_result: failed post request to Ruse resource creation: status {} error {}".format(
                    response.status_code, self._get_error(response)))
        resource_id = response.json()['id']
        self.resource_id = resource_id

        # upload the file to the resource
        uri = "{}/resource/{}/upload".format(ruse_server(), resource_id)
        with open(orig_filename, 'rb' if self.binary_file else 'r') as fh:
            files = {'file': fh}
            response = requests.put(uri, files=files)
            if response.status_code != requests.codes.ok and response.status_code != requests.codes.created:
                raise RuseServerException(
                    "upload_result: failed put request to Ruse resource id {} upload: status {} error {}".format(
                        response.status_code, resource_id, self._get_error(response)))

        return resource_id

    def resource_information(self, resource_id: int = None) -> None:
        """
        Retrieves information about the resource file from the ruse server

        :param resource_id: the resource id.  I None (default), the class attribute for resource_id should be set
        """

        if not resource_id:
            resource_id = self.resource_id
        assert resource_id

        uri = "{}/resource/{}".format(ruse_server(), resource_id)
        response = requests.get(uri)
        if response.status_code != requests.codes.ok:
            raise RuseServerException(
                "upload_result: failed get request to Ruse resource id {} upload: status {} error {}".format(
                    response.status_code, resource_id, self._get_error(response)))

        data = response.json()

        assert data['id'] == resource_id
        self.resource_id = resource_id
        for field in ['owner_id', 'name', 'orig_filename', 'uuid', 'ext']:
            setattr(self, field, data[field])

    def delete(self, resource_id: Optional[int]) -> None:
        if not resource_id:
            resource_id = self.resource_id
        assert resource_id
        uri = "{}/resource/{}".format(ruse_server(), resource_id)
        requests.delete(uri)

    def download(self) -> str:
        """
        Downloads the resource to the service runner directory: normally jsvcrunner will do this, but the manual
        download here is useful for testing

        :return: path to downloaded file
        """

        assert self.resource_id
        assert self.ext
        assert self.uuid

        work_dir = ruse_work_dir()
        file = os.path.abspath(os.path.join(work_dir, 'resources', "{}.{}".format(self.uuid, self.ext)))
        uri = "{}/resource/{}/download".format(ruse_server(), self.resource_id)
        response = requests.get(uri)
        with open(file, 'wb' if self.binary_file else 'w') as fh:
            fh.write(response.content)

        self.downloaded_file = file
        return file

    @classmethod
    def resource_map(cls, resources: Iterable) -> Dict[str, str]:
        """
        Create expected resource map json.  For service testing

        :param resources:
        :return: dictionary of resource ids mapped to files
        """

        return {str(r.resource_id): r.downloaded_file for r in resources}
