Attribute VB_Name = "DISQUE"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << DISQUE >>
'  Requis   :
'  Externes :
'
'  Descrip. : Disque Lecture/Ecriture de données
'
Option Explicit
Sub LireVec1(ByRef NomFichier As String, _
             ByRef n As Integer, _
             ByRef T() As Double)

    n = 0
    Open NomFichier For Input As #1
        Do While Not EOF(1)
            n = n + 1
            ReDim Preserve T(n)
            Input #1, T(n)
        Loop
    Close #1
End Sub
Sub LireVec2(ByRef NomFichier As String, _
             ByRef n As Integer, _
             ByRef T1() As Double, _
             ByRef T2() As Double)

    n = 0
    Open NomFichier For Input As #1
        Do While Not EOF(1)
            n = n + 1
            ReDim Preserve T1(n), T2(n)
            Input #1, T1(n), T2(n)
        Loop
    Close #1
End Sub
Sub LireVec3(ByRef NomFichier As String, _
             ByRef n As Integer, _
             ByRef T1() As Double, _
             ByRef T2() As Double, _
             ByRef T3() As Double)

    n = 0
    Open NomFichier For Input As #1
        Do While Not EOF(1)
            n = n + 1
            ReDim Preserve T1(n), T2(n), T3(n)
            Input #1, T1(n), T2(n), T3(n)
        Loop
    Close #1
End Sub
Sub LireVec4(ByRef NomFichier As String, _
             ByRef n As Integer, _
             ByRef T1() As Double, _
             ByRef T2() As Double, _
             ByRef T3() As Double, _
             ByRef T4() As Double)

    n = 0
    Open NomFichier For Input As #1
        Do While Not EOF(1)
            n = n + 1
            ReDim Preserve T1(n), T2(n), T3(n), T4(n)
            Input #1, T1(n), T2(n), T3(n), T4(n)
        Loop
    Close #1
End Sub
Sub LireMatNN(ByRef NomFichier As String, _
              ByRef n As Integer, _
              ByRef A() As Double, _
     Optional ByVal s As String = " ")

    Dim i As Integer, j As Integer, c As String
    Dim L As Integer, k As Integer, p As Integer
    
    s = Mid(Trim(s), 1, 1)
    If s = "" Then s = " "
    Open NomFichier For Input As #1
        Input #1, n
        ReDim A(n, n)
        For i = 1 To n
            Line Input #1, c
            For j = 1 To n - 1
                p = prv_Pos(c, s, A(i, j))
                If p = 0 Then Stop
            Next j
            A(i, j) = Val(c)
        Next i
    Close #1
End Sub
Sub LireMatNm(ByRef NomFichier As String, _
              ByRef n As Integer, _
              ByRef m As Integer, _
              ByRef A() As Double, _
     Optional ByVal s As String = " ")

    Dim i As Integer, j As Integer, c As String
    Dim L As Integer, k As Integer, p As Integer
    
    s = Mid(Trim(s), 1, 1)
    If s = "" Then s = " "
    Open NomFichier For Input As #1
        Input #1, n, m
        ReDim A(n, m)
        For i = 1 To n
            Line Input #1, c
            For j = 1 To m - 1
                p = prv_Pos(c, s, A(i, j))
                If p = 0 Then Stop
            Next j
            A(i, j) = Val(c)
        Next i
    Close #1
End Sub
Sub LireSys(ByRef NomFichier As String, _
            ByRef n As Integer, _
            ByRef A() As Double, _
            ByRef b() As Double, _
   Optional ByVal s As String = " ")

    Dim i As Integer, j As Integer, c As String
    Dim p As Integer
    
    s = Mid(Trim(s), 1, 1)
    If s = "" Then s = " "
    Open NomFichier For Input As #1
        Input #1, n
        ReDim A(n, n), b(n)
        For i = 1 To n
            Line Input #1, c
            For j = 1 To n
                p = prv_Pos(c, s, A(i, j))
                If p = 0 Then Stop
            Next j
            b(i) = Val(c)
        Next i
    Close #1
End Sub
Sub EcrireVec1(ByRef NomFichier As String, _
               ByVal n As Integer, _
               ByRef T() As Double)
    
    Dim i As Integer
    
    Open NomFichier For Output As #1
        For i = 1 To n
            Write #1, T(i)
        Next i
    Close #1
End Sub
Sub EcrireVec2(ByRef NomFichier As String, _
               ByVal n As Integer, _
               ByRef T1() As Double, _
               ByRef T2() As Double)

    Dim i As Integer
    
    Open NomFichier For Output As #1
        For i = 1 To n
            Write #1, T1(i), T2(i)
        Next i
    Close #1
End Sub
Sub EcrireVec3(ByRef NomFichier As String, _
               ByVal n As Integer, _
               ByRef T1() As Double, _
               ByRef T2() As Double, _
               ByRef T3() As Double)

    Dim i As Integer
    
    Open NomFichier For Output As #1
        For i = 1 To n
            Write #1, T1(i), T2(i), T3(i)
        Next i
    Close #1
End Sub
Sub EcrireVec4(ByRef NomFichier As String, _
               ByVal n As Integer, _
               ByRef T1() As Double, _
               ByRef T2() As Double, _
               ByRef T3() As Double, _
               ByRef T4() As Double)

    Dim i As Integer
    
    Open NomFichier For Output As #1
        For i = 1 To n
            Write #1, T1(i), T2(i), T3(i), T4(i)
        Next i
    Close #1
End Sub
Sub EcrireMatNN(ByRef NomFichier As String, _
                ByVal n As Integer, _
                ByRef A() As Double)

    Dim i As Integer, j As Integer, c As String
    Dim L As Integer, k As Integer, p As Integer
    
    Open NomFichier For Output As #1
        Write #1, n
        For i = 1 To n
            c = ""
            For j = 1 To n
                c = c + Str(A(i, j)) + "   "
            Next j
           Write #1, c
        Next i
    Close #1
End Sub
Sub EcrireMatNm(ByRef NomFichier As String, _
                ByVal n As Integer, _
                ByVal m As Integer, _
                ByRef A() As Double)

    Dim i As Integer, j As Integer, c As String
    Dim L As Integer, k As Integer, p As Integer
    
    Open NomFichier For Output As #1
        Write #1, n, m
        For i = 1 To n
            c = ""
            For j = 1 To m
                c = c + Str(A(i, j)) + "   "
            Next j
            Write #1, c
        Next i
    Close #1
End Sub
Sub EcrireSys(ByRef NomFichier As String, _
              ByVal n As Integer, _
              ByRef A() As Double, _
              ByRef b() As Double)

    Dim i As Integer, j As Integer, c As String
    Dim p As Integer
    
    Open NomFichier For Output As #1
        Write #1, n
        For i = 1 To n
            c = ""
            For j = 1 To n
                c = c + Str(A(i, j)) + "   "
            Next j
            c = c + Str(b(i))
            Write #1, c
        Next i
    Close #1
End Sub
Private Function prv_Pos(ByRef c As String, _
                         ByVal s As String, _
                         ByRef x As Double) As Integer

    Dim p As Integer, k As Integer, L As Integer
    
    c = Trim(c): L = Len(c)
    If L = 0 Then prv_Pos = 0: Exit Function
    p = 0
    k = 0
    Do
        k = k + 1
        If Mid(c, k, 1) = s Then p = k
    Loop Until k = L Or p <> 0
    If p <> 0 Then
        x = Val(Left(c, p - 1))
        c = Mid(c, p + 1, L - p)
    Else
        x = 0
    End If
    prv_Pos = p
End Function
