# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:31:08 2024

@author: dibas
"""

# Author: Marco Simoes
# Adapted from Java's implementation of Rui Pedro Paiva
# Teoria da Informacao, LEI, 2022

import sys
import numpy as np
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
        saida = []
        BFINAL = 0
        while not BFINAL == 1:

            BFINAL = self.readBits(1)

            BTYPE = self.readBits(2)
            if BTYPE != 2:
                print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks + 1))
                return

            # --- STUDENTS --- ADD CODE HERE
            
            # SEMANA 1
            print("\n")
            infos = self.getInfos()
            HLIT = infos[0]
            HDIST = infos[1]
            HCLEN = infos[2]
            print(HLIT)
            print(HDIST)
            print(HCLEN)
            print("\n")
            
            # Tabela para a ordem específica dos códigos
            code_length_order = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]

            # Array para armazenar comprimentos dos códigos (HCLEN + 4 valores, cada um com 3 bits)
            code_lengths = [0] * len(code_length_order)

            for i in range(HCLEN + 4):
                length = self.readBits(3)  # Lê 3 bits para o comprimento do código
                code_lengths[code_length_order[i]] = length  # Armazena na ordem especificada

            print("Array de comprimentos dos códigos:", code_lengths)
            print("\n")
            
            # SEMANA 2
            MAX_COMP = max(code_lengths)  # Pega o maior valor do comprimento de codigos
            array_cont_comp = self.contagemComprimentos(code_lengths, MAX_COMP)  # armazena em array_cont_comp as ocorrencias de cada comprimento
            print("Array soma da ocorrencia de cada comprimento: ", array_cont_comp)
            print("\n")
            
            arrayInicioCodigo = self.ArrayCodigos(array_cont_comp, MAX_COMP)  # armazena em arrayInicioCodigo o incio de codigo para cada comprimento
            print("Array para inicio de cada codigo: ", arrayInicioCodigo)
            print("\n")
            
            arrayCodigos = self.gerarCodigos(array_cont_comp, arrayInicioCodigo, MAX_COMP)
            print("Array dos codigos: ", arrayCodigos)
            print("\n")
            
            arrayIndices = self.gerarArrayIndices(code_lengths, MAX_COMP)
            print("Array de Índices:", arrayIndices)
            print("\n")
            
            hft = HuffmanTree()
            verbose = False
            for i, indice in enumerate(arrayIndices):
                hft.addNode(arrayCodigos[i], indice, verbose)
                
            # SEMANA 3
            array_lit_comp = self.comprimentoCodigos(HLIT, 257, hft)
            #print(array_lit_comp)
            print("\n")
            
            array_dist = self.comprimentoCodigos(HDIST, 1, hft)
            #print(array_dist)
            print("\n")
            
            #SEMANA 4
            
            # Codigos de huffman para os literais / comprimentos (repeticao das funçoes ja criadas)
            MAX_COMP = max(array_lit_comp)  # Pega o maior valor do comprimento de codigos
            array_cont_comp = self.contagemComprimentos(array_lit_comp, MAX_COMP)  # armazena em array_cont_comp as ocorrencias de cada comprimento
            print("Array soma da ocorrencia de cada comprimento: ", array_cont_comp)
            print("\n")
            
            arrayInicioCodigo = self.ArrayCodigos(array_cont_comp, MAX_COMP)  # armazena em arrayInicioCodigo o incio de codigo para cada comprimento
            print("Array para inicio de cada codigo: ", arrayInicioCodigo)
            print("\n")
            
            arrayCodigosLIT = self.gerarCodigos(array_cont_comp, arrayInicioCodigo, MAX_COMP)
            print("Array dos codigos: ", arrayCodigosLIT)
            print("\n")
            
            arrayIndicesLIT = self.gerarArrayIndices(array_lit_comp, MAX_COMP)
            print("Array de Índices:", arrayIndices)
            print("\n")
            
            # Codigos de huffman para as distancias (repeticao das funçoes ja criadas)
            MAX_COMP = max(array_dist)  # Pega o maior valor do comprimento de codigos
            array_cont_comp = self.contagemComprimentos(array_dist, MAX_COMP)  # armazena em array_cont_comp as ocorrencias de cada comprimento
            print("Array soma da ocorrencia de cada comprimento: ", array_cont_comp)
            print("\n")
            
            arrayInicioCodigo = self.ArrayCodigos(array_cont_comp, MAX_COMP)  # armazena em arrayInicioCodigo o incio de codigo para cada comprimento
            print("Array para inicio de cada codigo: ", arrayInicioCodigo)
            print("\n")
            
            arrayCodigosDIST = self.gerarCodigos(array_cont_comp, arrayInicioCodigo, MAX_COMP)
            print("Array dos codigos: ", arrayCodigosDIST)
            print("\n")
            
            arrayIndicesDIST = self.gerarArrayIndices(array_dist, MAX_COMP)
            print("Array de Índices:", arrayIndices)
            print("\n")
            
            # Arvore para os literais / comprimentos
            CLC = HuffmanTree()
            for i, indice in enumerate(arrayIndicesLIT):
                CLC.addNode(arrayCodigosLIT[i], indice, verbose)
            
            # Arvore para as distancias
            D = HuffmanTree()
            for i, indice in enumerate(arrayIndicesDIST):
                D.addNode(arrayCodigosDIST[i], indice, verbose)    
                
            saida = self.descompactacao(CLC, D, saida)
            
            #

            # update number of blocks read
            numBlocks += 1

        # SEMANA 5
        
        with open(self.gzh.fName, 'wb') as arquivo:
            arquivo.write(bytearray(saida))
        arquivo.close()

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
    
    # Gera os simbolos (indices) relacionados a cada codigo para a arvore de huffman
    def gerarArrayIndices(self, code_lengths, maxComp):
        indices = []
        for comprimento in range(1, maxComp + 1):  # Para cada comprimento válido
            for i, length in enumerate(code_lengths):
                if length == comprimento:
                    indices.append(i)  # Adiciona o índice do código com esse comprimento
        return indices
    
    def comprimentoCodigos(self, tamanho, tamanhoAdicional, hft):
        tamanhoTotal = tamanho + tamanhoAdicional
        arrayComprimentos = [0] * tamanhoTotal
        i = 0
        while i < tamanhoTotal:
            hft.resetCurNode()
            pos = -2
            
            while pos == -2:
                noAtual = self.readBits(1)
                pos = hft.nextNode(str(noAtual))
                
            if pos == -1:
                print("Error")
                return -1
                
            elif pos == 16:
                bitsAdicionais = self.readBits(2)
                numExtensaoValor = 3 + bitsAdicionais
                valorAnterior = arrayComprimentos[i - 1]
                fim = i + numExtensaoValor
                arrayComprimentos[i : fim] = [valorAnterior] * numExtensaoValor
                i = fim
                
            elif pos == 17:
                bitsAdicionais = self.readBits(3)
                numZerosAdicionais = 3 + bitsAdicionais
                fim = i + numZerosAdicionais
                arrayComprimentos[i : fim] = [0] * numZerosAdicionais
                i = fim
                
            elif pos == 18:
                bitsAdicionais = self.readBits(7)
                numZerosAdicionais = 11 + bitsAdicionais
                fim = i + numZerosAdicionais
                arrayComprimentos[i : fim] = [0] * numZerosAdicionais
                i = fim
                
            else:
                arrayComprimentos[i] = pos;
                i += 1
                
        return arrayComprimentos
            
    def descompactacao(self, CLC, DIST, saida):
        indice = 0
        
        # Arrays para a logica dos comprimentos
        arrayCodeComp = np.arange(257, 286)
        arrayBitsExtraComp = [0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0]
        arrayBaseComp = [3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258]
        
        # Arrays para a logica das distancias
        arrayCodeDist = np.arange(0, 30)
        arrayBitsExtraDist = [0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13]
        arrayBaseDist = [1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577]
        
        while True:
            CLC.resetCurNode()
            pos = -2
            
            while pos == -2:
                bit = self.readBits(1)
                pos = CLC.nextNode(str(bit))
            
            if pos == -1:
                print("Error")
                break
            
            elif pos < 256:
                literal = pos
                saida.append(literal)
                indice += 1
                
            elif pos == 256:
                print("Fim do bloco")
                break
            
            else:
                comp = self.decodifica_comp(pos, arrayCodeComp, arrayBaseComp, arrayBitsExtraComp)
                dist = self.decodifica_dist(DIST, arrayCodeDist, arrayBaseDist, arrayBitsExtraDist)
                
                for i in range(comp):
                    saida.append(saida[indice - dist + i])
                indice += comp
            
        return saida
    
    def decodifica_comp(self, pos, arrayCode, arrayBaseComp, arrayBitsEtras):
        for indice, elemento in enumerate(arrayCode):
            if elemento == pos:
                if pos >= 0 and pos <= 3:
                    comp = arrayBaseComp[indice]
                else:
                    bitsExtras = self.readBits(arrayBitsEtras[indice])
                    comp = arrayBaseComp[indice] + bitsExtras
        return comp
                
    def decodifica_dist(self, DIST, arrayCode, arrayBaseDist, arrayBitsEtras):
        DIST.resetCurNode()
        pos = -2
        
        while pos == -2:
            bit = self.readBits(1)
            pos = DIST.nextNode(str(bit))
        
        if pos == -1:
            print("Error")
            return
        
        for indice, elemento in enumerate(arrayCode):
            if pos == elemento:
                if pos >= 0 and pos <= 3:
                    distancia = arrayBaseDist[indice]
                else:
                    bitsExtras = self.readBits(arrayBitsEtras[indice])
                    distancia = arrayBaseDist[indice] + bitsExtras
        return distancia

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
        

if __name__ == "__main__":
    # Lista de arquivos a serem descompactados
    arquivos_gzip = ["FAQ.txt.gz", "sample_image.jpeg.gz", "sample_audio.mp3.gz", "sample_large_text.txt.gz"]
    
    # Processar cada arquivo da lista
    for arquivo in arquivos_gzip:
        try:
            print(f"Descompactando arquivo: {arquivo}")
            gzip_obj = GZIP(arquivo)  # Criar objeto GZIP
            gzip_obj.decompress()  # Chamar método para descompactar
            print(f"Arquivo {arquivo} descompactado com sucesso.\n")
        except Exception as e:
            print(f"Erro ao descompactar o arquivo {arquivo}: {e}\n")
