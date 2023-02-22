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
Sub Lecture1(ByRef NomFichier As String, _
             ByRef N As Integer, _
             ByRef VecX() As Double)

    N = 0
    Open NomFichier For Input As #1
        Do While Not EOF(1)
            N = N + 1
            ReDim Preserve VecX(N)
            Input #1, VecX(N)
        Loop
    Close #1
End Sub
Sub Lecture2(ByRef NomFichier As String, _
             ByRef N As Integer, _
             ByRef VecX() As Double, _
             ByRef VecY() As Double)

    N = 0
    Open NomFichier For Input As #1
        Do While Not EOF(1)
            N = N + 1
            ReDim Preserve VecX(N), VecY(N)
            Input #1, VecX(N), VecY(N)
        Loop
    Close #1
End Sub
Sub Ecriture1(ByRef NomFichier As String, _
              ByVal N As Integer, _
              ByRef VecX() As Double)
    
    Dim i As Integer
    
    Open NomFichier For Output As #1
        For i = 1 To N
            Write #1, VecX(i)
        Next i
    Close #1
End Sub
Sub Ecriture2(ByRef NomFichier As String, _
              ByVal N As Integer, _
              ByRef VecX() As Double, _
              ByRef VecY() As Double)

    Dim i As Integer
    
    Open NomFichier For Output As #1
        For i = 1 To N
            Write #1, VecX(i), VecY(i)
        Next i
    Close #1
End Sub

