# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:31:08 2024

@author: dibas
"""

# Author: Marco Simoes
# Adapted from Java's implementation of Rui Pedro Paiva
# Teoria da Informacao, LEI, 2022

import sys
from huffmantree import HuffmanTree


class GZIPHeader:
    ''' class for reading and storing GZIP header fields '''

    ID1 = ID2 = CM = FLG = XFL = OS = 0
    MTIME = []
    lenMTIME = 4
    mTime = 0

    # bits 0, 1, 2, 3 and 4, respectively (remaining 3 bits: reserved)
    FLG_FTEXT = FLG_FHCRC = FLG_FEXTRA = FLG_FNAME = FLG_FCOMMENT = 0

    # FLG_FTEXT --> ignored (usually 0)
    # if FLG_FEXTRA == 1
    XLEN, extraField = [], []
    lenXLEN = 2

    # if FLG_FNAME == 1
    fName = ''  # ends when a byte with value 0 is read

    # if FLG_FCOMMENT == 1
    fComment = ''  # ends when a byte with value 0 is read

    # if FLG_HCRC == 1
    HCRC = []

    def read(self, f):
        ''' reads and processes the Huffman header from file. Returns 0 if no error, -1 otherwise '''

        # ID 1 and 2: fixed values
        self.ID1 = f.read(1)[0]
        if self.ID1 != 0x1f: return -1  # error in the header

        self.ID2 = f.read(1)[0]
        if self.ID2 != 0x8b: return -1  # error in the header

        # CM - Compression Method: must be the value 8 for deflate
        self.CM = f.read(1)[0]
        if self.CM != 0x08: return -1  # error in the header

        # Flags
        self.FLG = f.read(1)[0]

        # MTIME
        self.MTIME = [0] * self.lenMTIME
        self.mTime = 0
        for i in range(self.lenMTIME):
            self.MTIME[i] = f.read(1)[0]
            self.mTime += self.MTIME[i] << (8 * i)

        # XFL (not processed...)
        self.XFL = f.read(1)[0]

        # OS (not processed...)
        self.OS = f.read(1)[0]

        # --- Check Flags
        self.FLG_FTEXT = self.FLG & 0x01
        self.FLG_FHCRC = (self.FLG & 0x02) >> 1
        self.FLG_FEXTRA = (self.FLG & 0x04) >> 2
        self.FLG_FNAME = (self.FLG & 0x08) >> 3
        self.FLG_FCOMMENT = (self.FLG & 0x10) >> 4

        # FLG_EXTRA
        if self.FLG_FEXTRA == 1:
            # read 2 bytes XLEN + XLEN bytes de extra field
            # 1st byte: LSB, 2nd: MSB
            self.XLEN = [0] * self.lenXLEN
            self.XLEN[0] = f.read(1)[0]
            self.XLEN[1] = f.read(1)[0]
            self.xlen = self.XLEN[1] << 8 + self.XLEN[0]

            # read extraField and ignore its values
            self.extraField = f.read(self.xlen)

        def read_str_until_0(f):
            s = ''
            while True:
                c = f.read(1)[0]
                if c == 0:
                    return s
                s += chr(c)

        # FLG_FNAME
        if self.FLG_FNAME == 1:
            self.fName = read_str_until_0(f)

        # FLG_FCOMMENT
        if self.FLG_FCOMMENT == 1:
            self.fComment = read_str_until_0(f)

        # FLG_FHCRC (not processed...)
        if self.FLG_FHCRC == 1:
            self.HCRC = f.read(2)

        return 0


class GZIP:
    ''' class for GZIP decompressing file (if compressed with deflate) '''

    gzh = None
    gzFile = ''
    fileSize = origFileSize = -1
    numBlocks = 0
    f = None

    bits_buffer = 0
    available_bits = 0

    def __init__(self, filename):
        self.gzFile = filename
        self.f = open(filename, 'rb')
        self.f.seek(0, 2)
        self.fileSize = self.f.tell()
        self.f.seek(0)

    def decompress(self):
        ''' main function for decompressing the gzip file with deflate algorithm '''

        numBlocks = 0

        # get original file size: size of file before compression
        origFileSize = self.getOrigFileSize()
        print(origFileSize)

        # read GZIP header
        error = self.getHeader()
        if error != 0:
            print('Formato invalido!')
            return

        # show filename read from GZIP header
        print(self.gzh.fName)

        # MAIN LOOP - decode block by block
        BFINAL = 0
        while not BFINAL == 1:

            BFINAL = self.readBits(1)

            BTYPE = self.readBits(2)
            if BTYPE != 2:
                print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks + 1))
                return

            # --- STUDENTS --- ADD CODE HERE
            
            # SEMANA 1
            infos = self.getInfos()
            HLIT = infos[0]
            HDIST = infos[1]
            HCLEN = infos[2]
            print(HLIT)
            print(HDIST)
            print(HCLEN)
            
            # Tabela para a ordem específica dos códigos
            code_length_order = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]

            # Array para armazenar comprimentos dos códigos (HCLEN + 4 valores, cada um com 3 bits)
            code_lengths = [0] * len(code_length_order)

            for i in range(HCLEN + 4):
                length = self.readBits(3)  # Lê 3 bits para o comprimento do código
                code_lengths[code_length_order[i]] = length  # Armazena na ordem especificada

            print("Array de comprimentos dos códigos:", code_lengths)
            
            # SEMANA 2
            MAX_COMP = max(code_lengths)  # Pega o maior valor do comprimento de codigos
            array_cont_comp = self.contagemComprimentos(code_lengths, MAX_COMP)  # armazena em array_cont_comp as ocorrencias de cada comprimento
            print("Array soma da ocorrencia de cada comprimento: ", array_cont_comp)
            
            arrayInicioCodigo = self.ArrayCodigos(array_cont_comp, MAX_COMP)  # armazena em arrayInicioCodigo o incio de codigo para cada comprimento
            print("Array para inicio de cada codigo: ", arrayInicioCodigo)
            
            arrayCodigos = self.gerarCodigos(array_cont_comp, arrayInicioCodigo, MAX_COMP)
            print("Array dos codigos: ", arrayCodigos)
            #

            # update number of blocks read
            numBlocks += 1

        # close file

        self.f.close()
        print("End: %d block(s) analyzed." % numBlocks)
        
    def getInfos(self):
        return [self.readBits(5), self.readBits(5), self.readBits(4)]
    
    def contagemComprimentos(self, comprimentos, maxComp):
        array = [0] * (maxComp + 1)  # Cria um array com tamanho do maior comprimento + 1, para que inclua tambem esse elemento no array
        for i in range(maxComp + 1):
            count = 0
            for j in comprimentos:
                if(j == i):
                    count += 1
            array[i] = count
        return array
            
    def ArrayCodigos(self, contagens, maxComp):
        code = 0
        contagens[0] = 0
        next_code = [0] * (maxComp + 1)
        for bits in range(1, maxComp + 1):
            code = (code + contagens[bits - 1]) << 1 # Logica para o começo do proximo codigo de comprimento (equivalente a multiplica po 2)
            next_code[bits] = code;
        return next_code
    
    def gerarCodigos(self, arrayContagens, arrayInicio, maxComp):
        soma = sum(arrayContagens)
        arrayCodigos = [0] * soma  # Criação do array de códigos, com o tamanho das ocorrências (excluindo o 0)
        indice = 0
        for i in range(1, maxComp + 1):
            if arrayContagens[i] != 0:
                inicioCodigo = arrayInicio[i]  # Começo do código de comprimento i
                fimCodigo = arrayContagens[i] + inicioCodigo  # Fim do código de comprimento i
                for j in range(inicioCodigo, fimCodigo):
                    code = bin(j)[2:]  # Remove o prefixo "0b"
                
                    # Completa o código com zeros à esquerda para atingir o comprimento desejado
                    if len(code) < i:
                        code = code.zfill(i)  # preenche com zeros à esquerda até atingir o tamanho `i`
                    arrayCodigos[indice] = code
                    indice += 1
        return arrayCodigos

    def getOrigFileSize(self):
        ''' reads file size of original file (before compression) - ISIZE '''

        # saves current position of file pointer
        fp = self.f.tell()

        # jumps to end-4 position
        self.f.seek(self.fileSize - 4)

        # reads the last 4 bytes (LITTLE ENDIAN)
        sz = 0
        for i in range(4):
            sz += self.f.read(1)[0] << (8 * i)

        # restores file pointer to its original position
        self.f.seek(fp)

        return sz

    def getHeader(self):
        ''' reads GZIP header'''

        self.gzh = GZIPHeader()
        header_error = self.gzh.read(self.f)
        return header_error

    def readBits(self, n, keep=False):
        ''' reads n bits from bits_buffer. if keep = True, leaves bits in the buffer for future accesses '''

        while n > self.available_bits:
            self.bits_buffer = self.f.read(1)[0] << self.available_bits | self.bits_buffer
            self.available_bits += 8

        mask = (2 ** n) - 1
        value = self.bits_buffer & mask

        if not keep:
            self.bits_buffer >>= n
            self.available_bits -= n

        return value
        

if __name__ == '__main__':

    # gets filename from command line if provided
    fileName = "FAQ.txt.gz"
    if len(sys.argv) > 1:
        fileName = sys.argv[1]

    # decompress file
    gz = GZIP(fileName)
    gz.decompress()
